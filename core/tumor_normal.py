import collections
import os
import subprocess
import string
import warnings
import scipy.stats
import numpy as np


class Variant(object):
    ''' Store information about a variant
    '''
    def __init__(self,*args,**kwargs):
        '''
        '''
        self.refUMI = None
        self.altUMI = None
        self.oddsRatio = None
        self.pval = None
        self.adjPval = None
        

def parseVariant(row,vmt_col):
    ''' Helper function to parse multi-allelic variants and their VMT
    from cut/anno.vcf/txt files

    :param list row
    :return ((chrom,pos,ref,alt),vmt)
    '''
    if row[4].find(",")!=-1: # multi-allelic
        alts = row[4].split(",")
        vmts = row[vmt_col].split(",")
        for i,alt in enumerate(alts):
            yield ((row[0],row[1],row[3],alt),int(vmts[i])) # ((chrom,pos,ref,alt),vmt)
    else:
        yield ((row[0],row[1],row[3],row[4]),int(row[vmt_col]))

        
def parseAnnoTxt(readSet,refUMICounts):
    ''' Parse annotated txt flat file , storing variant information
    :param str  readSet
    :param dict refUMICounts
    '''
    variantObj = collections.defaultdict(Variant)
    with open(readSet + ".smCounter.anno.txt", "r") as IN:
        for line in IN:
            contents = line.strip("\n").split("\t")
            if contents[0] == "CHROM": # header"
                continue
            
            for variantKey,altUMICount in parseVariant(contents,11):
                pos    = variantKey[1]
                variantObj[variantKey].refUMI = refUMICounts[pos]
                variantObj[variantKey].altUMI = altUMICount
                
    return variantObj

    
def parseSmCounterAllFile(readSet):
    ''' Parse smCounter all file , storing ref  UMI count for each base
    Enforce pval cutoffs if needed
    :param str  readSet
    :param bool enforceCutoffs
    '''
    cutoff = 6
    minCutoff = {"INDEL":2,"SNP":2} ## Cutoff for low PI file

    ## dicts for storing base position -> refUMI
    allVars = {}

    # parse smCounter all.txt file
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

            allVars[POS] = refUMI
                
    return allVars


def fet(tumorVars, normalVars):
    ''' Perform fisher's exact test
    :param dict tumorVars :
    :param dict normalVars :
    '''
    for var in tumorVars:
        refUMITumor  =  tumorVars[var].refUMI
        altUMITumor  =  tumorVars[var].altUMI
        refUMINormal =  normalVars[var].refUMI
        altUMINormal =  normalVars[var].altUMI
        
        # do F.E.T
        oddsRatio,pval = scipy.stats.fisher_exact([[refUMITumor, refUMINormal], [altUMITumor, altUMINormal]])
        tumorVars[var].pval      = pval
        tumorVars[var].oddsRatio = oddsRatio

    return (tumorVars)


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


def applyTNFilter(f,ftype,tumorVarsFiltered,pValCutoff):
    ''' Add information from the F.E.T test to output files
    :param str  f
    :param str  ftype
    :param dict tumorVarsFiltered
    :param float pValCutoff
    '''
    # conditional to ignore header
    headerIgnore = {
        "cutVcf"  : lambda line:line.startswith("#"),
        "cutTxt"  : lambda line:line.startswith("CHROM"),
        "annoVcf"  : lambda line:line.startswith("#"),
        "annoTxt"  : lambda line:line.startswith("CHROM")        
    }
    assert ftype in headerIgnore, "tumor_normal: Invalid file type : {}".format(ftype)
    
    # filter column
    filterCol = {
        "cutVcf"  :  6,
        "cutTxt"  : -1,
        "annoVcf" :  6,
        "annoTxt" :  6
    }

    def simpleParseVariant(row,ftype):
        ''' 
        Helper function to parse multi-allelic variants from cut/anno.vcf/txt files
        Similar to parseVariant, except does not return vmt
        :param list row
        :param ftype i,e, cutTxt,cutVcf,annoTxt,annoVcf
        :return (chrom,pos,ref,alt)
        '''
        if ftype == 'cutTxt':
            refCol = 2            
            altCol = 3
        else:
            refCol = 3
            altCol = 4            
            
        if row[altCol].find(",")!=-1: # multi-allelic
            alts = row[altCol].split(",")
            for i,alt in enumerate(alts):
                yield (row[0],row[1],row[refCol],alt) # (chrom,pos,ref,alt)
        else:
            yield (row[0],row[1],row[refCol],row[altCol])


    f = f[:-1] if f.endswith('/') else f # remove trailing /
    tempFile = f + ".temp"
    with open(f,"r") as IN, open(tempFile,"w") as OUT:
        for line in IN:
            contents = line.strip("\n").split("\t")
            if headerIgnore[ftype](line):
                OUT.write(line)
                continue
            filterField = contents[filterCol[ftype]]
            tnFilter = []
            for variant in simpleParseVariant(contents,ftype):
                if variant in tumorVarsFiltered: # had enough UMIs for F.E.T and was present in normal sample
                    pval      = tumorVarsFiltered[variant].pval
                    adjPval   = tumorVarsFiltered[variant].adjPval
                    oddsRatio = tumorVarsFiltered[variant].oddsRatio
                    assert pval is not None, "tumor_normal: Logical Error for variant : {}".format(variant)
                    assert adjPval is not None, "tumor_normal: Logical Error for variant : {}".format(variant)
                    assert oddsRatio is not None, "tumor_normal: Logical Error for variant : {}".format(variant)                    
                    if adjPval < pValCutoff:
                        tnFilter.append("TN/oddsRatio:%e/pval:%e/pvalAdj:%e"%(oddsRatio,pval,adjPval))
            # update filter
            if tnFilter:
                filterField = filterField + ";" + ",".join(tnFilter)

            contents[filterCol[ftype]] = filterField
            OUT.write("\t".join(contents))
            OUT.write("\n")

    # replace the input file with the updated new one
    subprocess.check_call("mv {temp} {out}".format(temp=tempFile,out=f),shell=True)
    

def tumorNormalVarFilter(cfg):
    ''' Filter Tumor variants
    '''
    readSetNormal =  cfg.readSetMatchedNormal
    readSetTumor  =  cfg.readSet
    umiCutoff     =  int(cfg.umiCutoff) # umi cutoff for F.E.T
    pValCutoff    =  float(cfg.pValCutoff) # cutoff for adjusted p values from F.E.T
    print("\ntumor_normal: Started filtering Tumor variants")

    # do nothing if zero variants from tumor read set
    if not os.path.isfile(readSetTumor + ".smCounter.cut.txt"):
        return

    # parse all files for refUMI count
    normalRefUMICounts = parseSmCounterAllFile(readSetNormal)
    tumorRefUMICounts  = parseSmCounterAllFile(readSetTumor)
    
    # parse tumor anno txt, add refUMICounts to obj
    tumorVars  = parseAnnoTxt(readSetTumor,tumorRefUMICounts)
    # parse normal anno txt, add refUMICounts to obj
    normalVars = parseAnnoTxt(readSetNormal,normalRefUMICounts)
    
    # filter variants with < 10 ref + alt UMI
    tumorVarsTmp  =  {k:v for k,v in tumorVars.items() if v.refUMI + v.altUMI >= 10}
    normalVarsTmp =  {k:v for k,v in normalVars.items() if v.refUMI + v.altUMI >=10}

    # filter tumor variants not in normal
    tumorVarsFiltered  = {k:v for k,v in tumorVars.items() if k in normalVarsTmp}
    normalVarsFiltered = {k:v for k,v in normalVars.items() if k in tumorVarsFiltered} # do not need to do this, done for consistency, we won't use normalVariants for any outputs

    # perform F.E.T and store pvalue and oddsRatio
    tumorVarsFiltered = fet(tumorVarsFiltered,normalVarsFiltered)    

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

    # parse and update anno,cut.vcf/txt
    files = [readSetTumor + ".smCounter.cut.vcf", readSetTumor + ".smCounter.cut.txt",
             readSetTumor + ".smCounter.anno.vcf", readSetTumor + ".smCounter.anno.txt"]
    fileTypes = ["cutVcf","cutTxt","annoVcf","annoTxt"]
    for i,f in enumerate(files):
        applyTNFilter(f, fileTypes[i], tumorVarsFiltered,pValCutoff)
    print("tumor_normal: Finished filtering Tumor variants")
    
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
