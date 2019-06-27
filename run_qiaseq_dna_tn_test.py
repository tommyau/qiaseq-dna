import ConfigParser
import sys
import multiprocessing
# our modules
import core.run_log
import core.run_config
#import core.prep
import core.align
#import core.umi_filter
#import core.umi_mark
#import core.umi_merge
import core.primer_clip
import core.samtools
import core.tumor_normal
#import core.sm_counter_wrapper
#import metrics.sum_specificity
#import metrics.sum_uniformity_primer
#import metrics.sum_primer_umis
#import metrics.sum_all
#import metrics.umi_frags
#import metrics.umi_depths
#import misc.process_ion
#import misc.tvc
#import annotate.vcf_complex
#import annotate.vcf_annotate

#--------------------------------------------------------------------------------------
# call input molecules, build consenus reads, align to genome, trim primer region
#--------------------------------------------------------------------------------------
def run(args):
    readSet, paramFile, vc = args
 
    # read run configuration file to memory
    cfg = core.run_config.run(readSet,paramFile)   
   
 
def run_tumor_normal(readSet,paramFile,vc):
    ''' Wrapper around run() for tumor-normal analysis
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

    # initialize logger
    #core.run_log.init(normal)    
    #run((normal,paramFile,vc))
    # close log file
    #core.run_log.close()

    # initialize logger
    core.run_log.init(tumor)
    #run((tumor,paramFile,vc))    
    ## Additional analysis steps
    cfg = core.run_config.run(tumor,paramFile)
    import os
    import shutil
    # restore from backups if they exist
    if os.path.exists(tumor + ".smCounter.cut.txt.bak"):
        shutil.copyfile(tumor + ".smCounter.cut.txt.bak",
                        tumor + ".smCounter.cut.txt")
        shutil.copyfile(normal + ".smCounter.cut.txt.bak",
                        normal + ".smCounter.cut.txt")
        shutil.copyfile(tumor + ".smCounter.all.txt.bak",
                        tumor + ".smCounter.all.txt")
        shutil.copyfile(normal + ".smCounter.all.txt.bak",
                        normal + ".smCounter.all.txt")
        
    core.tumor_normal.tumorNormalVarFilter(cfg, normal, tumor)
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
        
        # initialize logger
        core.run_log.init(readSet)
        
        run((readSet,paramFile,vc))
        cfg = core.run_config.run(readSet,paramFile)
        core.tumor_normal.runCopyNumberEstimates(cfg)
        
        # close log file
        core.run_log.close() 
