import collections
import os
import subprocess
import string

import scipy
import numpy as np

def parseSmCounterAllFile(readSet,enforceCutoff):
    ''' Parse smCounter all file , storing ref and alt UMI count for each variant
    Enforce cutoffs for tumor variants, return type changes based on this param

    :param str  readSet
    :param bool enforceCutoff
    '''
    cutoff = 6
    minCutoff = {"INDEL":2,"SNP":2} ## Cutoff for low PI file

    # init namedtuple obj
    fields = ["refUMI", "altUMI", "oddsRatio", "pval", "adjPval"]
    variantInfo = collections.namedtuple('variantInfo', fields)
    # initialize field values to None by default
    variantInfo.__new__.__defaults__ = (None,) * len(variantInfo._fields)
    
    ## some dicts for storing variants
    # passing log(pval) threshold
    tumorAboveCutoffVars = collections.defaultdict(variantInfo)
    # passing lowQ log(pval) threshold
    tumorLowQVars = collections.defaultdict(variantInfo)
    # all variants
    allVars = collections.defaultdict(variantInfo)

    seen = set() # keep track of variants encountered for debug
    # parse smCounter all.txt file, keeping lowQ and above threshold variants if tumor sample, otherwise keep everything
    with open(readSet + ".smCounter.all.txt", "r") as IN:
        for line in IN:
            contents = line.strip("\n").split("\t")
            
            if contents[0] == "CHROM": # header
                continue

            CHROM, POS, REF, ALT, TYPE, sUMT, sForUMT, sRevUMT, sVMT, sForVMT, sRevVMT, sVMF, sForVMF, sRevVMF, VDP, VAF, RefForPrimer, RefRevPrimer, primerOR, pLowQ, hqUmiEff, allUmiEff, refMeanRpb, altMeanRpb, rpbEffectSize, repType, hpInfo, simpleRepeatInfo, tandemRepeatInfo, DP, FR, MT, UFR, sUMT_A, sUMT_T, sUMT_G, sUMT_C, logpval, FILTER = line.strip().split("\t") # using same variable notation as smCounter vcf module for ease of understanding

            if TYPE == "0":
                continue
            
            if ALT == "DEL":
                continue

            QUAL = logpval if logpval != "NA" else 0.00

            try:
                fQUAL = float(QUAL)
            except ValueError:
                fQUAL = 0.00

            if enforceCutoffs: # used for tumor sample
                if fQUAL < minCutoff[TYPE.upper()]:
                    continue
            
            key = tuple(contents[0],contents[1],contents[3],contents[4]) # chrom,pos,ref,alt
            assert key not in seen, "tumor_normal: Duplicate Variant Encountered !"
            
            if REF[0] == "A": # there could be deletions, pysam pileup uses first base as reference
                refUMI = sUMT_A
            elif REF[0] == "T":
                refUMI = sUMT_T
            elif REF[0] == "G":
                refUMI = sUMT_G
            elif REF[0] == "C":
                refUMI = sUMT_C
            else:
                raise Exception("tumor_normal: LogicalError. Could not discern correct reference base and umi count.")
            
            if enforceCutoffs: # used for tumor sample
                if fQUAL < cutoff: # lowQ
                    lowQVars[key].refUMI = refUMI
                    lowQVars[key].altUMI = sVMT
                else:
                    aboveCutoffVars[key].refUMI = refUMI
                    aboveCutoffVars[key].altUMI = sVMT
            else:
                allVars[key].refUMI = refUMI
                allVars[key].altUMI = sVMT
                
    if enforceCutoffs:
        return (aboveCutoffVars, lowQVars)
    else:
        return allVars
    
def fet(tumorAboveCutoffVars, tumorLowQVars, normalVars, umiCutoff):
    ''' Perform fisher's exact test
    :param dict tumorAboveCutoffVars :
    :param dict tumorLowQVars :
    :param dict normalVars :
    :param int  umiCutoff :
    '''
    for var in tumorAboveCutoffVars:
        refUMITumor  =  tumorAboveCutoffVars[var].refUMI
        altUMITumor  =  tumorAboveCutoffVars[var].altUMI
        refUMINormal =  normalVars[var].refUMI
        altUMINormal =  normalVars[var].altUMI
        
        # do F.E.T only if UMI count > cutoff
        if refUMITumor + altUMITumor > umiCutoff and refUMINormal + altUMINormal > umiCutoff:
            pval, oddsRatio = scipy.stats.fisher_exact([[refUMITumor, refUMINormal], [altUMITumor, altUMINormal]])
            tumorAboveCutoffVars[var].pval      = pval
            tumorAboveCutoffVars[var].oddsRatio = oddsRatio

    for var in tumorLowQVars:
        refUMITumor, altUMITumor    =  tumorLowQVars[var]
        refUMINormal, altUMINormal  =  normalVars[var]

        # do F.E.T only if UMI count > cutoff
        if refUMITumor + altUMITumor > umiCutoff and refUMINormal + altUMINormal > umiCutoff:
            pval, oddsRatio = scipy.stats.fisher_exact([[refUMITumor, refUMINormal], [altUMITumor, altUMINormal]])
            tumorLowQVars[var].pval      = pval
            tumorLowQVars[var].oddsRatio = oddsRatio
            
    return (tumorAboveCutoffVars, tumorLowQVars)

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
    
def tumorNormalVarFilter(cfg):
    ''' Filter Tumor variants
    '''
    readSetNormal =  cfg.readSetMatchedNormal
    readSetTumor  =  cfg.readSet
    umiCutoff     =  cfg.umiCutoff # umi cutoff for F.E.T
    
    # do nothing if zero variants from tumor read set
    if not os.path.isfile(readSetTumor + ".smCounter.cut.txt"):
        return

    # parse tumor file
    tumorAboveCutoffVars, tumorLowQVars = parseSmCounterAllFile(readSetTumor,enforceCutoffs=True)
    alltumorVars = set(tumorAboveCutoffVars.keys() + tumorLowQVars.keys())
    # parse normal file
    normalVars = parseSmCounterAllFile(readSetTumor)

    # filter variants with < 10 ref + alt UMI
    tumorAboveCutoffVarsFiltered  =  {k:v for k,v in tumorAboveCutoffVars if v.refUMI + v.altUMI >= 10}
    tumorLowQVarsFiltered         =  {k:v for k,v in tumorLowQVars if v.refUMI + v.altUMI >= 10}
    normalVarsFiltered            =  {k:v for k,v in normalVars if v.refUMI + v.altUMI >=10}

    # perform F.E.T and store pvalue and oddsRatio
    tumorAboveCutoffVarsFiltered, tumorLowQVarsFiltered = fet(tumorAboveCutoffVarsFiltered,
                                                              tumorLowQVarsFiltered,
                                                              normalVarsFiltered)

    allTumorVarsFiltered   =  tumorAboveCutoffVarsFiltered.keys() + tumorLowQVarsFiltered.keys()
    allTumorPvals  =  [v.pval for k,v in tumorAboveCutoffVarsFiltered] + \
                      [v.pval for k,v in tumorlowQVarsFiltered]
    # correct for FDR using Benjamini-Hotchberg
    adjPvals = benjaminiHotchbergCorrection(allTumorPvals)
    for i,var in enumerate(allTumorVars):
        if var in tumorAboveCutoffVarsFiltered:
            tumorAboveCutoffVarsFiltered[var].adjPval = adjPvals[i]
        elif var in tumorLowQVarsVarsFiltered:
            tumorLowQVarsVarsFiltered[var].adjPval = adjPvals[i]
        else:
            raise Exception("tumor_normal: Logical Error. Mimsatch in variants after FDR pval adjustment")

    # parse and update cut.vcf/txt and belowThreshold lowQ file
    f = readSetTumor + ".smCounter.cut.vcf"
    ftype = "cutVcf"
    applyTNFilter(f, tumorAboveCutoffVarsFiltered, allTumorVars, ftype)
    
    f = readSetTumor + ".smCounter.cut.txt"
    ftype = "cutTxt"
    applyTNFilter(f, tumorAboveCutoffVarsFiltered, allTumorVars, ftype)
    
    f = readSetTumor + ".smCounter.lowQ.txt"
    ftype = "lowQ"
    applyTNFilter(f, tumorLowQVarsFiltered, allTumorVars, ftype)
    
def applyTNFilter(relevantVars,allTumorVars,f,ftype):
    ''' Add information from the F.E.T test to output files
    :param str  f
    :param dict relevantVars
    :param dict allTumorVars
    :param str  ftype
    '''
    # conditional to ignore header
    headerIgnore = {
        "cutVcf"  : lambda line:line.startswith("#"),
        "cutTxt"  : lambda line:line.startswith("CHROM"),
        "lowQ"    : lambda line:line.find("CHROM")!=-1
    }
    assert ftype in headerIgnore, "tumor_normal: Invalid file type : {}".format(ftype)
    
    # filter column
    filterCol = {
        "cutVcf"  : 6,
        "cutTxt"  : 6,
        "lowQ"    : -1
    }
    # return variant key tuple from row
    def parseVariant(ftype):
        ''' Parse row based on file type
        '''
        if ftype == "lowQ":
            yield (contents[1],contents[2],contents[3],contents[4])
        else:
            if contents[4].find(",")!=-1: # multi-allelic
                alt1,alt2 = contents[4].split(",")
                yield (contents[0],contents[1],contents[3],alt1)
                yield (contents[0],contents[1],contents[3],alt2)

    f = f[:-1] if f.endswith('/') else f # remove trailing /
    tempFile = f + ".temp"
    with open(f,"r") as IN, open(tempFile,"w") as OUT:
        for line in IN:
            contents = line.strip("\n").split("\t")
            if headerIgnore[ftype](line):
                OUT.write(line)
                continue
            filter = contents[filterCol(ftype)]
            tnFilter = []
            for variant in parseVariant(contents):
                assert variant in allTumorVars, "tumor_normal: Logical Error. Variant from all file not present in : {} file".format(ftype)
                if variant in relevantVars: # had enough UMIs for F.E.T
                    pval    = round(relevantVars.pval,3)
                    adjPval = round(relevantVars.adjPval,3)
                    tnFilter.append("-".join(["TN","pval:{p}","qval:{q}"]))
            # update filter
            if fetFilter:
                filter = filter + ";" + ",".join(tnFilter)
            contents[filterCol(ftype)] = filter
            OUT.write("\t".join(contents))
            OUT.write("\n")

    # replace the input file with the updated new one
    subprocess.check_call("mv {temp} {in}".format(temp=tempFile,out=f),shell=True)
    
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
