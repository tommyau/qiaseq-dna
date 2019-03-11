import ConfigParser
import sys
import multiprocessing
import shutil
# our modules
import core.run_log
import core.run_config
import core.prep
import core.align
import core.umi_filter
import core.umi_mark
import core.umi_merge
import core.primer_clip
import core.samtools
import core.tumor_normal
import core.sm_counter_wrapper
import metrics.sum_specificity
import metrics.sum_uniformity_primer
import metrics.sum_primer_umis
import metrics.sum_all
import metrics.umi_frags
import metrics.umi_depths
import misc.process_ion
import misc.tvc
import annotate.vcf_complex
import annotate.vcf_annotate

#--------------------------------------------------------------------------------------
# call input molecules, build consenus reads, align to genome, trim primer region
#--------------------------------------------------------------------------------------
def run(args,tumorNormal):
    readSet, paramFile, vc = args
    # initialize logger
    if not tumorNormal:
        core.run_log.init(readSet)
 
    # read run configuration file to memory
    cfg = core.run_config.run(readSet,paramFile)

    # trim adapters , umi and primers (this module spawns multiple processes)
    core.prep.run(cfg)

    readFileIn1 = readSet + ".prep.R1.fastq"
    readFileIn2 = readSet + ".prep.R2.fastq"
    bamFileOut  = readSet + ".align.bam"    
    if cfg.platform.lower() == "illumina":
        # align trimmed reads to genome using BWA MEM
        core.align.run(cfg, readFileIn1, readFileIn2, bamFileOut)
    else: # use tmap for ion torrent reads        
        misc.process_ion.alignToGenomeIon(cfg, readFileIn1, bamFileOut)
  
    # call putative unique input molecules using BOTH UMI seq AND genome alignment position on random fragmentation side
    
    bamFileIn  = readSet + ".align.bam"
     
    core.umi_filter.run(cfg, bamFileIn)
    core.umi_mark.run(cfg)   
    metrics.umi_frags.run(cfg)
       
    metrics.umi_depths.run(cfg,vc)   
    core.umi_merge.run(cfg, bamFileIn)
    
    # soft clip primer regions from read alignments
    bamFileIn  = readSet + ".umi_merge.bam"
    bamFileOut = readSet + ".primer_clip.bam"
    core.primer_clip.run(cfg,bamFileIn,bamFileOut,False)
 
    # additional metrics to generate
    metrics.sum_primer_umis.run(cfg) # primer-level umi and read metrics
    metrics.sum_specificity.run(cfg) # priming specificity
    metrics.sum_uniformity_primer.run(cfg) # primer-level uniformity

    # sort the final BAM file, to prepare for downstream variant calling
    bamFileIn  = readSet + ".primer_clip.bam"
    bamFileOut = readSet + ".bam"
    core.samtools.sort(cfg,bamFileIn,bamFileOut)   
   
    if cfg.duplex.lower() == "false": # do not run smCounter for duplex reads
 
        if cfg.platform.lower() != "illumina": # ion reads
            misc.tvc.run(cfg)
   
        # run smCounter variant calling
        numVariants = core.sm_counter_wrapper.run(cfg, paramFile, vc)
        
        if cfg.platform.lower() != "illumina":
            numVariants = misc.tvc.smCounterFilter(cfg,vc)
   
        # create complex variants, and annotate using snpEff
        if not tumorNormal:
            post_smcounter_work(numVariants, readSet, cfg)
            # close log file
            core.run_log.close()
 
def post_smcounter_work(numVariants, readSet, cfg):
    ''' Additional Steps after smCounter
    :param int numVariants
    :param str readSet
    :param lambda obj cfg
    '''
    if numVariants > 0:
        # convert nearby primitive variants to complex variants
        bamFileIn  = readSet + ".bam"
        vcfFileIn  = readSet + ".smCounter.cut.vcf"
        vcfFileOut = readSet + ".smCounter.cplx.vcf"
        annotate.vcf_complex.run(cfg, bamFileIn, vcfFileIn, vcfFileOut, vc)
            
        # annotate variants in the VCF file
        vcfFileIn  = readSet + ".smCounter.cplx.vcf"
        vcfFileOut = readSet + ".smCounter.anno.vcf"
        annotate.vcf_annotate.run(cfg, vcfFileIn, vcfFileOut,vc)

    else: # create a header only anno.vcf from cut.vcf 
        vcfFileIn  = readSet + ".smCounter.cut.vcf"
        vcfFileOut = readSet + ".smCounter.anno.vcf"
        shutil.copyfile(vcfFileIn,vcfFileOut)
        
        # aggregate all metrics
        metrics.sum_all(cfg)

def run_tumor_normal(readSet,paramFile,vc):
    ''' Wrapper for tumor-normal analysis
    '''
    # 2 read set names which are space delimited
    readSets = filter(None,readSet.split(" "))
    assert len(readSets) == 2, "Tumor-Normal Analysis requires exactly 2 read sets !"
 
    # read parameter file
    parser = ConfigParser.SafeConfigParser()
    parser.optionxform = str
    parser.read(paramFile)
 
    tumor = None
    normal = None
    for section in parser.sections():
        if section not in ['general','smCounter']:
            for (paramName, paramVal) in parser.items(section):
                if paramName == 'sampleType' and paramVal.lower() == 'normal':
                    normal = section
                elif paramName == 'sampleType' and paramVal.lower() == 'tumor':
                    tumor = section
     
    assert tumor!=None and normal!=None, "Could not sync read set names supplied with config file !"
    core.run_log.init(tumor)
    run((tumor,paramFile,vc),tumorNormal=True)
    print("--"*20)
    run((normal,paramFile,vc),tumorNormal=True)

    ## Compare Tumor Normal variants and update filter
    if vc == 'v2': # use new TN filter
        print("--"*20)
        cfg = core.run_config.run(tumor,paramFile)
        core.tumor_normal.tumorNormalVarFilter(cfg)
        print("--"*20)

    ## Create cplx,anno.txt/vcf and sum.all files
    numVariants = core.tumor_normal.getNumVariants(tumor)
    post_smcounter_work(numVariants,tumor, cfg)
    print("--"*20)
    numVariants = core.tumor_normal.getNumVariants(normal)
    cfg = core.run_config.run(normal,paramFile)
    post_smcounter_work(numVariants,normal, cfg)
    print("--"*20)

    ## Run old variant substraction code if using v1
    if vc == 'v1':
        print("--"*20)
        print('Warning: Doing naive substraction of normal variants from tumor. Please use smCounter-v2  for newer Tumor-Normal variant Filter')
        cfg = core.run_config.run(tumor,paramFile)
        core.tumor_normal.removeNormalVariants(cfg)

    ## Run Quandico for CNV    
    core.tumor_normal.runCopyNumberEstimates(cfg)

    # close log file
    core.run_log.close()

#-------------------------------------------------------------------------------------
# main program for running from shell 
#-------------------------------------------------------------------------------------
if __name__ == "__main__":

    if len(sys.argv) > 6 :
        print "\nRun as : python run_qiaseq_dna.py <param_file> <v1/v2> <single/tumor-normal> <readSet(s)>\n"
        sys.exit(-1)
  
    paramFile = sys.argv[1]
    vc = sys.argv[2]
    analysis = sys.argv[3]
    readSet   = " ".join(sys.argv[4:]) # 2 readSets in case of tumor-normal
 
    if analysis.lower() == "tumor-normal":      
        run_tumor_normal(readSet,paramFile,vc)
    else: # Single sample, might still need to run quandico
        run((readSet,paramFile,vc),tumorNormal = False)
        cfg = core.run_config.run(readSet,paramFile)
        core.tumor_normal.runCopyNumberEstimates(cfg)
