import collections
import copy
import operator
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
smcounter_v2 = __import__("qiaseq-smcounter-v2.vcf")

# some constants, learned from existing datasets
_a_ = 0.381774 # het vmf lower bound
_b_ = 0.597383 # het vmf upper bound
_c_ = 0.949296 # homozygous threshold when pure tumor, e.g. 0.99


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

        # for HomRef Tumor Variants, store the corresponding normal variant key
        self.normalVarKey = None

    @staticmethod
    def isHet(vmf):
        ''' Returns if a variant is Heterozygous
            Used for identifying normal het sites

        :param float vmf: The variant allele frequency
        :rtype: bool
        '''
        return _a_ <= vmf <= _b_

    @staticmethod
    def isHomAlt(tumorVmf, normalVmf, p):
        ''' Returns if a variant is Homozygous for the alt allele
            Used for tumor variants

        :param float tumorVmf: The tumor variant allele frequency
        :param float normalVmf: The corresponding normal variant allele frequency
        :param float p: Tumor Purity
        :rtype: bool
        '''
        return tumorVmf > (p * _c_ + (1 - p) * normalVmf)

    @staticmethod
    def isHomRef(tumorVmf, normalVmf, p):
        ''' Returns if a variant is tending to Homozygous for the ref allele
            Used for tumor variants
        :param float tumorVmf: The tumor variant allele frequency
        :param float normalVmf: The corresponding normal variant allele frequency
        :param float p: Tumor Purity
        :rtype: bool
        '''
        return tumorVmf < (p * (1 - _c_) + (1 - p) * normalVmf)

    def classifyLOH(self, nvmf, p):
        ''' Check if the tumor variant is because of Loss of Heterozygocity
        :param float nvmf: The normal variant allele frequency
        :param float p: Tumor Purity
        '''
        if self.isHomAlt(self.vmf, nvmf, p) and self.isHet(nvmf):
            self.varclass = 'LOH'
            assert self.vmf > nvmf, "Logical Error !"

        elif self.isHomRef(self.vmf, nvmf, p) and self.isHet(nvmf):
            assert self.vmf < nvmf, "Logical Error !"
            self.varclass = 'LOH_HomRef'


    def classifyGermlineOrSomatic(self, nvmf, cutoff):
        ''' Classify a tumor variant as Germline of somatic
        :param float nvmf: The normal variant allele frequency
        :param float cutoff : adjusted pValue cutoff from FET
        '''
        if self.adjPval >= cutoff:
            self.varclass = 'Germline_Risk'

        else:
            assert self.vmf > nvmf, "Logical Error !"
            self.varclass = 'Somatic'


def getNumVariants(readSet):
    ''' Helper function to return number of variants in cut.txt file
    :param str readSet: The readSet name
    '''
    n = 0
    with open(readSet + '.smCounter.cut.txt', 'r') as IN:
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
    allVars    = collections.defaultdict(Variant) # key : (chrom, pos, ref, alt)
    homRefVars = collections.defaultdict(list)    # key : (chrom, pos, ref) 
    
    # parse smCounter all.txt file
    with open(readSet + ".smCounter.all.txt", "r") as IN:
        for line in IN:
            contents = line.strip("\n").split("\t")

            if contents[0] == "CHROM": # header
                continue

            CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = contents # using same variable notation as smCounter vcf module for ease of understanding

            key = (CHROM, POS, REF, ALT)

            if sVMF == '.': # ignore sites with no UMI coverage
                continue

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
            elif len(REF) == 1 and REF[0] == "N":
                continue # ignore sites with reference base N
            else:
                raise Exception("tumor_normal: LogicalError. Could not discern correct reference base and umi count.")
            allVars[key].refUMI  = refUMI
            allVars[key].altUMI  = int(sVMT)
            allVars[key].vmf     = sVMF
            allVars[key].qual    = fQUAL
            allVars[key].vartype = TYPE.upper()

            sUMT = int(sUMT)
            refAlleleFreq = float(refUMI)/sUMT
            if(refAlleleFreq >= _c_ and fQUAL < 2.00): # check if variant is homozygous ref - To capture sites with < 3 Alt UMIs
                key2 = (CHROM, POS, REF)
                homRefVars[key2].append(copy.deepcopy(allVars[key]))

    return (allVars, homRefVars)


def handleTumorHomRefVars(normalAllVars, tumorAllVars, tumorHomRefVars):
    ''' Create a dummy variant key in tumorAllVars for :
    Homozygous Reference Tumor sites where the Normal Sample has a Heterozygous Variant.
    The alt allele is different in the Tumor and Normal at these sites
    Store a dummy variant at this site with the same alt allele as the Normal
    '''
    found = False

    for variantKey in normalAllVars:

        if variantKey not in tumorAllVars:

            if Variant.isHet(normalAllVars[variantKey].vmf): # check if normal is Het
                chrom, pos, ref, alt = variantKey
                key2 = (chrom, pos, ref)

                if key2 in tumorHomRefVars: # check if tumor is Hom Ref
                    found = True
                    
                    # if > 2 variants at this site, pick the one with highest variant UMI count
                    sortedAlts = sorted([(variant.altUMI,variant) for variant in \
                                         tumorHomRefVars[key2]],
                                        key = operator.itemgetter(0),
                                        reverse = True)

                    variant = sortedAlts[0][1]
                    tumorAllVars[variantKey] = variant
                    for i in range(len(tumorHomRefVars[key2])):
                        tumorHomRefVars[key2][i].normalVarKey = variantKey
                    # if normal has 2 het sites , could get two possible variants here to store, current logic will keep only last one, but the LOH filter at this site should still be captured in the tumor.
                    # an example : Normal : A->C (38%vmf) ; A->T (38%vmf)
                    #              Tumor  : A->G (1%vmf, i.e. 99% ref, QUAL should be < 2 as well)


def fet(variantKey, tumorAllVars, normalAllVars):
    ''' Perform fisher's exact test
    :param tuple variantKey, i.e. (chrom,pos,ref,alt)
    '''
    if variantKey not in tumorAllVars:
        return (-1, -1, variantKey)

    # log(pval) cutoff
    cutoff = 6

    qualTumor    =  tumorAllVars[variantKey].qual

    if qualTumor >= cutoff:

        refUMITumor  =  tumorAllVars[variantKey].refUMI
        altUMITumor  =  tumorAllVars[variantKey].altUMI
        refUMINormal =  normalAllVars[variantKey].refUMI
        altUMINormal =  normalAllVars[variantKey].altUMI

        # do F.E.T
        oddsRatio, pval = scipy.stats.fisher_exact(
            [[refUMITumor, refUMINormal], [altUMITumor, altUMINormal]])
        return (pval, oddsRatio, variantKey)
    else:
        return (-1, -1, variantKey)


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


def updateFilter(readSet, outFile, tumorVarsFiltered, tumorHomRefVars, tumorAllVars, isTumor):
    ''' Update Filter column for variants in smCounter all.txt file
    Add a new column with adjPval info
    :param str readSet
    :param str outFile
    :param defaultdict(Variant Obj) tumorVarsFiltered : All variants where TN FET test was done
    :param defaultdict(list) tumorHomRefVars : Homozygous Reference tumor variants
    :param defaultdict(Variant Obj) tumorAllVars : All parsed tumor variants 
    :param bool isTumor
    '''
    OUT = open(outFile, "w")
    OUT_debug = open(readSet + ".tumor_normal.all.txt", "w")

    for row in iterateSmCounterAllFile(readSet):
        if row[0] == "CHROM":
            row.append('TNFetPval')
            OUT.write("\t".join(row))
            OUT.write("\n")
            OUT_debug.write("\t".join(row))
            OUT_debug.write("\n")
            continue

        variantKey = (row[0], row[1], row[2], row[3])
        fltr       = row[-1]
        
        adjPval = '.'
        oddsRatio = '.'
        
        if isTumor and variantKey not in tumorVarsFiltered:
            key2 = (row[0], row[1], row[2])
            if key2 in tumorHomRefVars:
                normalKey = tumorHomRefVars[key2][0].normalVarKey # updated variant key
                if normalKey is None: # this site was not evaluated
                    #print("debug:variant : {} was not evaluated".format(variantKey))
                    pass
                else:
                    variantKey = normalKey # update variantKey to point to normal Het variant used to represent this site

        if variantKey in tumorVarsFiltered:
            adjPval   = tumorVarsFiltered[variantKey].adjPval
            oddsRatio = tumorVarsFiltered[variantKey].oddsRatio
            varclass  = tumorVarsFiltered[variantKey].varclass
            assert varclass is not None, "Unexpected Variant Filter!"
            assert (adjPval is not None and oddsRatio is not None) or varclass.startswith('LOH'), "Unexpected Variant info!"

            if isTumor:
                if varclass != 'Somatic':
                    newFltr = fltr + ';' + varclass if fltr != 'PASS' \
                              else varclass
                else:
                    newFltr = fltr

                row[-1] = newFltr
        else:
            assert tumorAllVars[variantKey].pval is None and tumorAllVars[variantKey].varclass is None, "LogicalError for variant: {}".format(variantKey)
        
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

        if variantKey in tumorVarsFiltered:
            OUT_debug.write("\t".join(row))
            OUT_debug.write("\n")

        OUT.write("\t".join(row))
        OUT.write("\n")

    OUT.close()
    OUT_debug.close()


def tumorNormalVarFilter(cfg, normal, tumor):
    ''' Filter Tumor variants
    '''
    readSetNormal =  normal
    readSetTumor  =  tumor
    umiCutoff     =  int(cfg.umiCutoff) # umi cutoff for F.E.T
    pValCutoff    =  float(cfg.pValCutoff) # cutoff for adjusted p values from F.E.T
    tumorPurity   =  float(cfg.tumorPurity)

    print("\ntumor_normal: Started filtering Tumor variants")

    # do nothing if zero variants from tumor read set
    if not os.path.isfile(readSetTumor + ".smCounter.cut.txt"):
        return

    # parse smCounter all.txt files
    normalAllVars, normalHomRefVars = parseSmCounterAllFile(readSetNormal)
    del(normalHomRefVars) # don't need these variants
    tumorAllVars, tumorHomRefVars = parseSmCounterAllFile(readSetTumor)
    # match up alt alleles with normal at tumor hom ref sites
    # normalAllVars and tumorAllVars will be updated
    handleTumorHomRefVars(normalAllVars, tumorAllVars, tumorHomRefVars)

    # flag LOH variants first
    LOHVars = set()
    for key in normalAllVars:
        # check if variant in tumor
        if key in tumorAllVars:
            nvmf = normalAllVars[key].vmf
            tumorAllVars[key].classifyLOH(nvmf, tumorPurity)
            if tumorAllVars[key].varclass is not None: # flagged as LOH
                LOHVars.add(key)

    # perform F.E.T and store pvalue and oddsRatio
    # This takes about 30 secs for ~ 1500 variants
    # I have left the fet function call as below to make it amenable to bring
    # in multiprocessing again if needed.
    # Need to declare normalAllVars, tumorAllVars and tumorHomRefVars as globals

    # I tried using multiprocessing.manager but could not get it to work
    # as class objects are not picklable, even namedtuple is not

    # Decided to leave it single threaded, ~ 30 secs speedup is not worth it
    # to have globals and clutter up code.
    
    for key in normalAllVars:
        if key not in LOHVars: # skip LOH variants
            pval, oddsRatio, key = fet(key, tumorAllVars, normalAllVars)
            if(pval != -1): # either tumor qual below cutoff, or variant not present in tumor all file
                tumorAllVars[key].pval = pval
                tumorAllVars[key].oddsRatio = oddsRatio

    # filter entries with no pvals
    tumorVarsFiltered  = {k:v for k,v in tumorAllVars.items() if v.pval is not None} # - variants flagged as LOH not here !

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
    for variantKey in tumorVarsFiltered:
        tumorVarsFiltered[variantKey].classifyVar(
            normalAllVars[variantKey].vmf, pValCutoff)

    # filter again - this time keep all variants which have a varclass, i.e. flagged as either LOH/LOH_HomRef/Germline_Risk/Somatic
    tumorVarsFiltered  = {k:v for k,v in tumorAllVars.items() if v.varclass is not None}

    # parse and update all.txt file
    tempFile1 = readSetTumor + ".smCounter.all.temp.txt"
    updateFilter(readSetTumor, tempFile1, tumorVarsFiltered, tumorHomRefVars, tumorAllVars, isTumor = True)
    tempFile2 = readSetNormal + ".smCounter.all.temp.txt"
    updateFilter(readSetNormal, tempFile2, tumorVarsFiltered, tumorHomRefVars, tumorAllVars, isTumor = False)

    # backup smCounter all files
    shutil.copyfile(readSetTumor + ".smCounter.cut.txt",
                    readSetTumor + ".smCounter.cut.txt.bak")
    shutil.copyfile(readSetNormal + ".smCounter.cut.txt",
                    readSetNormal + ".smCounter.cut.txt.bak")
    shutil.copyfile(readSetTumor + ".smCounter.all.txt",
                    readSetTumor + ".smCounter.all.txt.bak")
    shutil.copyfile(readSetNormal + ".smCounter.all.txt",
                    readSetNormal + ".smCounter.all.txt.bak")

    # re-run smCounter vcf creation module
    smcounter_v2.vcf.makeVcf(
        './', tempFile1, readSetTumor, cfg.genomeFile, isDuplex = False, tumorNormal = True)
    smcounter_v2.vcf.makeVcf(
        './', tempFile2, readSetNormal, cfg.genomeFile, isDuplex = False, tumorNormal = True)


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
    

def removeNormalVariants(cfg, normal, tumor):
    ''' Remove normal variants from tumor vcf
    '''
    readSetNormal = normal
    readSetTumor  = tumor

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
