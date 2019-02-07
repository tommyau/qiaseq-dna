import collections
import os
import math
import multiprocessing
import numpy as np
import shutil
import string
import scipy.stats
import subprocess
import sys
import warnings

# our modules
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
smcounter_v2 = __import__("qiaseq-smcounter-v2.vcf")

# some constants, learned from existing datasets
_a_ = 0.381774 # het vmf lower bound
_b_ = 0.597383 # het vmf upper bound
_c_ = 0.949296 # homozygous threshold when pure tumor, e.g. 0.99

# some globals used in function fet(), need this for multiprocessing
_normalAllVars_ = None
_tumorAllVars_  = None


class Variant(object):
    ''' Store information about a variant
    '''
    def __init__(self, *args, **kwargs):
        ''' Class Constructor
        '''
        self.refUMI      = None
        self.altUMI      = None
        self.qual        = None
        self.vmf         = None
        self.vartype     = None

        self.oddsRatio   = None
        self.pval        = None
        self.adjPval     = None
        self.varclass    = None

    def classifyVar(self, nvmf, p, cutoff):
        ''' Used for tumor variants
        :param float nvmf : Normal variant UMI allele frequency
        :param float p    : Tumor Purity
        :param float cutoff : pValue cutoff for FET

        :returns Variant classification, i.e. Germline_Risk, Somatic, LOH (Loss of Heterozygocity)
        :rtype string
        '''
        if self.adjPval >= cutoff:
            self.varclass = 'Germline_Risk'
        
        else:
            if self.vmf > nvmf:
                if (self.vmf > (p * _c_ + (1 - p) * nvmf) and \
                    (_a_ <= nvmf <= _b_)):
                    self.varclass = 'LOH'
                else:
                    self.varclass = 'Somatic'
            else:                
                if self.vmf < (p * (1 - _c_) + (1 - p) * nvmf) and \
                   _a_ <= nvmf <= _b_:
                    self.varclass = 'LOH_HomRef'
                else:
                    self.varclass = 'Germline_Risk'


def getNumVariants(readSet):
    ''' Helper function to return number of variants in cut.txt file
    :param str readSet: The readSet name
    '''
    n = 0
    with open(readSet + '.cut.txt', 'r') as IN:
        for line in IN:
            n+=1

    return n-1

def iterateSmCounterAllFile(readSet):
    ''' Yield line by line
    :param str readSet
    '''
    # parse smCounter all.txt file
    with open(readSet + ".smCounter.all.txt", "r") as IN:
        for line in IN:
            contents = line.strip("\n").split("\t")
            yield contents
            continue

def parseSmCounterAllFile(readSet):
    ''' Parse smCounter all file
    Enforce pval cutoffs if needed
    :param str  readSet
    '''
    # container to store variants
    allVars = collections.defaultdict(Variant)
    
    # parse smCounter all.txt file
    with open(readSet + ".smCounter.all.txt", "r") as IN:
        for line in IN:
            contents = line.strip("\n").split("\t")

            if contents[0] == "CHROM": # header
                continue

            CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = contents # using same variable notation as smCounter vcf module for ease of understanding

            key = (CHROM, POS, REF, ALT)

            sVMF = float(sVMF)/100

            QUAL = logpval if logpval != 'NA' else '0.00'

            try:
                fQUAL =  float(QUAL)
            except ValueError:
                fQUAL = 0.00

            if TYPE == "0":
                continue
            
            if ALT == "DEL":
                continue

            if REF[0] == "A": # there could be deletions, pysam pileup uses first base as reference
                refUMI = int(sUMT_A)
            elif REF[0] == "T":
                refUMI = int(sUMT_T)
            elif REF[0] == "G":
                refUMI = int(sUMT_G)
            elif REF[0] == "C":
                refUMI = int(sUMT_C)
            else:
                raise Exception("tumor_normal: LogicalError. Could not discern correct reference base and umi count.")
            allVars[key].refUMI  = refUMI
            allVars[key].altUMI  = sVMT
            allVars[key].vmf     = sVMF
            allVars[key].qual    = fQUAL
            allVars[key].vartype = TYPE.upper()


    return allVars


def fet(variant_key):
    ''' Perform fisher's exact test
    :param tuple variant_key, i.e. (chrom,pos,ref,alt)
    '''
    if variant_key not in _tumorAllVars_:
        return (-1, -1, variant_key)

    # log(pval) cutoffs
    cutoff = 6
    minCutoff = {"INDEL":2,"SNP":2} ## Cutoff for low PI file

    qualNormal   =  _normalAllVars_[variant_key].qual
    qualTumor    =  _tumorAllVars_[variant_key].qual

    if qualNormal >= minCutoff[_tumorAllVars_[variant_key].vartype] or \
       qualTumor >= minCutoff[_normalAllVars_[variant_key].vartype]:

        refUMITumor  =  _tumorAllVars_[variant_key].refUMI
        altUMITumor  =  _tumorAllVars_[variant_key].altUMI
        refUMINormal =  _normalAllVars_[variant_key].refUMI
        altUMINormal =  _normalAllVars_[variant_key].altUMI
   
        # do F.E.T
        oddsRatio, pval = scipy.stats.fisher_exact(
            [[refUMITumor, refUMINormal], [altUMITumor, altUMINormal]])
        return (pval, oddsRatio, variant_key)
    else:
        return (-1, -1, variant_key)


def benjaminiHotchbergCorrection(pvals):
    ''' Perform FDR correction using Benjamini-Hotchberg procedure
    :param list pvals
    :return : pvals corrected for FDR
    '''
    temp = np.array(pvals)
    sorted_indices = np.argsort(temp)
    m = len(sorted_indices) # number of tests
    rank = np.array(range(1,m+1))
    tempPAdjusted = m*temp[sorted_indices]/rank
    ## final adjustment
    ## we make sure that if a smaller q-value
    ## shows up in the list we set all of
    ## the q-values corresponding to smaller p-values
    ## to that corresponding q-value.
    
    ## We can do this by reversing the list
    ## Keeping a running tally of the minimum value
    ## and then re-reversing
    ## please see this helpful link : https://stackoverflow.com/questions/10323817/r-unexpected-results-from-p-adjust-fdr
    pAdjusted = np.minimum.accumulate(tempPAdjusted[::-1])[::-1]
    return pAdjusted[np.argsort(sorted_indices)]


def updateFilter(readSet, outFile, tumorVarsFiltered, isTumor):
    ''' Update Filter column for variants in smCounter all.txt file
    Add a new column with adjPval info
    :param str readSet
    :param str outFile
    :param Variant Obj tumorVarsFiltered
    :param bool isTumor
    '''
    OUT = open(outFile, 'w')
    
    for row in iterateSmCounterAllFile(readSet):
        if row[0] == "CHROM":
            row.append('TNFetPval')
            OUT.write("\t".join(row))
            OUT.write("\n")
            continue

        variantKey = (row[0], row[1], row[2], row[3])
        fltr      = row[-1]
        
        adjPval = '.'
        oddsRatio = '.'
        if variantKey in tumorVarsFiltered:
            adjPval   = tumorVarsFiltered[variantKey].adjPval
            oddsRatio = tumorVarsFiltered[variantKey].oddsRatio
            assert adjPval != None and oddsRatio != None, "Unexpected Variant FET info!"

            if isTumor:                
                varclass  = tumorVarsFiltered[variantKey].varclass
                assert varclass != None, "Unexpected Variant Filter!"

                if varclass != 'Somatic':
                    newFltr = fltr + ';' + varclass if fltr != 'PASS' \
                              else varclass
                else:
                    newFltr = fltr

                row[-1] = newFltr
        
        if adjPval != '.':
            # convert pvalue to phred scale
            adjPval   = max(1e-200, adjPval) # set min at 1e-200 to avoid log(0); same as smCounter
            adjPval   = round(-10 * math.log(adjPval,10), 2)
            adjPval   = abs(adjPval) # for weird -0.00

            oddsRatio = round(oddsRatio, 2)
            TNFetInfo = "%.2f;%.2f"%(adjPval,oddsRatio)
            TNFetPval = "%.2f"%(adjPval)
        else:
            TNFetInfo = "%s;%s"%(adjPval,oddsRatio)
            TNFetPval = "%s"%(adjPval)            
        
        row.append(TNFetPval)
        OUT.write("\t".join(row))
        OUT.write("\n")

    OUT.close()


def tumorNormalVarFilter(cfg):
    ''' Filter Tumor variants
    '''
    global _normalAllVars_
    global _tumorAllVars_

    readSetNormal =  cfg.readSetMatchedNormal
    readSetTumor  =  cfg.readSet
    umiCutoff     =  int(cfg.umiCutoff) # umi cutoff for F.E.T
    pValCutoff    =  float(cfg.pValCutoff) # cutoff for adjusted p values from F.E.T
    tumorPurity   =  float(cfg.tumorPurity)

    print("\ntumor_normal: Started filtering Tumor variants")

    # do nothing if zero variants from tumor read set
    if not os.path.isfile(readSetTumor + ".smCounter.cut.txt"):
        return

    # parse smCounter all.txt files
    _normalAllVars_ = parseSmCounterAllFile(readSetNormal)
    _tumorAllVars_  = parseSmCounterAllFile(readSetTumor)
    
    # perform F.E.T in parallel and store pvalue and oddsRatio
    pool = multiprocessing.Pool(int(cfg.numCores))
    for result in pool.imap_unordered(
            fet,_normalAllVars_.iterkeys(), chunksize=128):
        pval, oddsRatio, key = result
        if(pval != -1): # either tumor or normal var below cutoff, or variant not present in tumor all file
            _tumorAllVars_[key].pval = pval
            _tumorAllVars_[key].oddsRatio = oddsRatio
    # clear finished pool
    pool.close()
    pool.join()

    # filter entries with no pvals
    tumorVarsFiltered  = {k:v for k,v in _tumorAllVars_.items() if v.pval != None}

    # free some mem
    del(_tumorAllVars_)

    # store ordered pvals and variants
    tumorPvals   = []
    tumorVarKeys = []
    for key in tumorVarsFiltered:
        val = tumorVarsFiltered[key]
        tumorPvals.append(val.pval)
        tumorVarKeys.append(key)        

    # correct for FDR using Benjamini-Hotchberg
    adjPvals = benjaminiHotchbergCorrection(tumorPvals)    
    for i,var in enumerate(tumorVarKeys):
        tumorVarsFiltered[var].adjPval = adjPvals[i]

    # free some mem
    del(tumorPvals)
    del(tumorVarKeys)

    # update varClass
    for variant_key in tumorVarsFiltered:
        tumorVarsFiltered[variant_key].classifyVar(
            _normalAllVars_[variant_key].vmf, tumorPurity, pValCutoff)

    # parse and update all.txt file
    tempFile1 = readSetTumor + ".smCounter.all.temp.txt"
    updateFilter(readSetTumor, tempFile1, tumorVarsFiltered, isTumor = True)
    tempFile2 = readSetNormal + ".smCounter.all.temp.txt"
    updateFilter(readSetNormal, tempFile2, tumorVarsFiltered, isTumor = False)

    # backup smCounter all files
    shutil.copyfile(readSetTumor + ".smCounter.all.txt",
                    readSetTumor + ".smCounter.all.txt.bak")
    shutil.copyfile(readSetNormal + ".smCounter.all.txt",
                    readSetNormal + ".smCounter.all.txt.bak")

    # re-run smCounter vcf creation module
    smcounter_v2.vcf.makeVcf(
        './', tempFile1, readSetTumor, cfg.genomeFile, tumorNormal = True)
    smcounter_v2.vcf.makeVcf(
        './', tempFile2, readSetNormal, cfg.genomeFile, tumorNormal = True)


def runCopyNumberEstimates(cfg):
    ''' Run CNV analysis using quandico
    '''
    if cfg.runCNV.lower() == "false":
        return
    referenceUmiFiles = cfg.refUmiFiles.split(",")
    assert len(referenceUmiFiles) >= 1, "No reference UMI Files supplied !"

    # read reference UMI counts, normalize by read set, and get median for each primer across read sets (not necessary for only one reference read set)
    umiCountsAll = {}
    for fileName in referenceUmiFiles:
        # read MT counts from disk file
        umiCounts = []
        umiCountsTotal = 0
        for line in open(fileName,'r'):
            if not line.startswith("read set|"):
                vals = line.strip().split("|")
                (readSet, primer, strand, chrom, loc5, loc3, umiCount) = vals[0:7]
                (strand, loc5, umiCount) = map(int,(strand, loc5, umiCount))
                key = (chrom,strand,loc5,primer)
                umiCounts.append((key,umiCount))
                umiCountsTotal += umiCount

        # normalize to 1,000 mean MT depth and save in multi-readset hash
        meanUmiDepth = float(umiCountsTotal) / len(umiCounts)
        for key , umiCount in umiCounts:
            if key in umiCountsAll:
                vec = umiCountsAll[key]
            else:
                vec = []
                umiCountsAll[key] = vec
            vec.append(int(round(1000.00 * umiCount / meanUmiDepth)))

    # need to complement first base of negative strand primers
    dnaComplementTranslation = string.maketrans("ATGC", "TACG")

    # open output file for quandico input
    readSet = cfg.readSet
    umiFileReference = readSet + ".copy-number.reference.txt"
    fileout = open(umiFileReference, "w")

    # for each primer, get median MT depth across all reference read sets and write to disk
    for key, vec in umiCountsAll.iteritems():
        # get median MT depth for this primer
        idx = len(vec) / 2
        if len(vec) % 2 == 1: # odd length
            umiCount = vec[idx]
        else:
            umiCount = int(round((vec[idx-1] + vec[idx]) / 2.00))

        # unpack primer info
        (chrom, strand, loc5, primer) = key

        # write output in format needed by quandico
        refGenomeBase = primer[0]
        if strand == 1:
            refGenomeBase = refGenomeBase.translate(dnaComplementTranslation)
        outvec = (chrom, loc5+1, strand, refGenomeBase, "foobar", umiCount)
        fileout.write("\t".join((str(x) for x in outvec)))
        fileout.write("\n")

    # done writing input reference file
    fileout.close()

    # convert sample readset MT counts to format needed by quandico (NOTE: not normalized - quandico will do that)
    umiFileSample = readSet + ".copy-number.sample.txt"
    fileout = open(umiFileSample,"w")
    fileName = readSet + ".sum.primer.umis.txt"
    for line in open(fileName,"r"):
        if line.startswith("read set|"):
            continue
        vals = line.strip().split("|")
        (readSet_, primer, strand, chrom, loc5, loc3, mtCount) = vals[0:7]
        (strand, loc5, mtCount) = map(int,(strand, loc5, mtCount))
        refGenomeBase = primer[0]
        if strand == 1:
            refGenomeBase = refGenomeBase.translate(dnaComplementTranslation)
        outvec = (chrom, loc5+1, strand, refGenomeBase, "foobar", mtCount)
        fileout.write("\t".join((str(x) for x in outvec)))
        fileout.write("\n")
    fileout.close()

    # make work directory
    if not os.path.exists("_quandico_work_"):
        os.mkdir("_quandico_work_")

    # get code dir for quandico
    codeDirQuandico = cfg.quandicoDir

    # get reference genome path
    genomeFile = cfg.genomeFile + ".fai"

    # call quandico (note: hack in quandico.pl looks for ".copy-number" in the -b parameter)
    filePrefix = "{}.copy-number".format(readSet)
    cmd = "perl {}quandico.pl ".format(codeDirQuandico) \
        + "-E {}cluster.pl ".format(codeDirQuandico) \
        + "-y {}R/ ".format(codeDirQuandico) \
        + "-t _quandico_work_ " \
        + "-s data={} -s x=2 -s y=0 ".format(umiFileSample)  \
        + "-r data={} -r x=2 -r y=0 ".format(umiFileReference)  \
        + "-G data={} -G name=GRCh37 ".format(genomeFile) \
        + "-b {} ".format(filePrefix) \
        + " > {}.log 2>&1 ".format(filePrefix)
    print(cmd)
    subprocess.check_call(cmd, shell=True)
    

def removeNormalVariants(cfg):
    ''' Remove normal variants from tumor vcf
    '''
    readSetNormal = cfg.readSetMatchedNormal
    readSetTumor = cfg.readSet

    # do nothing if zero variants from tumor read set
    if not os.path.isfile(readSetTumor + ".smCounter.anno.txt"):
        return

    # save normal variants from flat detail file
    normalVariants = set()
    fileName = readSetNormal + ".smCounter.anno.txt"
    if os.path.isfile(fileName):
        with open(fileName,'r') as IN:
            for line in open(fileName,'r'):
                if not line.startswith("CHROM"):
                    vals = line.strip().split("\t")
                    normalVariants.add(tuple(vals[0:5]))

    # filter tumor VCF and corresponding tumor flat file
    for fileSuffix in (".smCounter.anno.txt",".smCounter.anno.vcf"):
        with open(readSetTumor + fileSuffix + ".temp",'w') as OUT:
            with open(readSetTumor + fileSuffix,'r') as IN:
                for line in IN:
                    if line.startswith("#"):
                        OUT.write(line)
                    else:
                        vals = line.strip().split('\t')
                        key = tuple(vals[0:5])
                        if key not in normalVariants:
                            OUT.write(line)
     # replace tumor vcf and flat file with the updated temp files
    for fileSuffix in (".smCounter.anno.txt",".smCounter.anno.vcf"):
        os.system("mv {temp} {f}".format(temp = readSetTumor + fileSuffix + ".temp",
                                         f    = readSetTumor + fileSuffix))
