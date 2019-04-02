import glob
import gzip
import os
import math
import subprocess
from multiprocessing.dummy import Pool as ThreadPool

# our modules
import primer_trim_ion

#-------------------------------------------------------------------------------------
def runShellCommand(cmd):
    # run shell command, capture stdout and stderror (assumes log is not large and not redirected to disk)
    try:
        log = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        error = False
    except subprocess.CalledProcessError, ex:
        log = ex.output
        error = True
        
    # print command log output
    for line in log.split("\n"):
        print("prep_trim: " + line.strip())
     
    # re-raise exception now that error detail has been printed
    if error:
        raise(ex)

#-------------------------------------------------------------------------------------
def worker(cmd):
    try:
        subprocess.check_call(cmd, shell=True)
        return True
    except:
        return False
  
#-------------------------------------------------------------------------------------
def splitReadFile(readFile,filePrefixOut,readSide,numBatchesMax,deleteLocalFiles):

    # validate existence of input read file
    if not os.path.isfile(readFile):
        raise UserWarning("input read file does not exist: " + readFile)
           
    # get number of reads in the fastq file
    isGzFile = readFile.endswith(".fastq.gz")
    if isGzFile:
        line = subprocess.check_output("zcat {} | wc -l".format(readFile), shell=True)
        numLines = int(line)
    else:
        line = subprocess.check_output("wc -l " + readFile, shell=True)
        vals = line.split(" ")
        numLines = int(vals[0])
    if numLines % 4 != 0 or numLines == 0:
        raise UserWarning("truncated or empty fastq read file {}".format(readFile))
    numReads = numLines / 4
    
    # set batch size
    batchSize = 4 * int(math.ceil(1.00 * numReads / numBatchesMax))
    
    # open input fastq read file
    if isGzFile:
        fileIn = gzip.open(readFile,"rb")
    else:
        fileIn = open(readFile,"r")
  
    # write fastq batches to disk - one batch to be done on each CPU core
    numLinesOut = 0
    batchNum = 0
    fileOut = None
    for line in fileIn:
        # open new file if needed
        if numLinesOut == 0:
            fileOut = open("{}.{:04d}.{}.fastq".format(filePrefixOut,batchNum,readSide), "w")
            
        # write fastq line to disk
        fileOut.write(line)
        numLinesOut += 1
           
        # close output file if done with batch
        if numLinesOut == batchSize:
            fileOut.close()
            numLinesOut = 0
            batchNum += 1
   
    # close out last batch
    if numLinesOut != 0:
        fileOut.close()
        batchNum += 1
    numBatches = batchNum
    
    # delete local input file if no longer needed
    if deleteLocalFiles and len(os.path.dirname(readFile)) == 0:
        os.remove(readFile)
     
    # done
    return numReads, numBatches

#-------------------------------------------------------------------------------------
def run(cfg):
    '''
    '''
    readSet          = cfg.readSet
    readFile1        = cfg.readFile1
    readFile2        = cfg.readFile2
    cutadaptDir      = cfg.cutadaptDir
    deleteLocalFiles = cfg.deleteLocalFiles
    numCores         = int(cfg.numCores)
    tagNameUmiSeq    = cfg.tagNameUmiSeq
    tagNamePrimer    = cfg.tagNamePrimer
    tagNamePrimerErr = cfg.tagNamePrimerErr
    primerFile       = cfg.primerFile
    primer3Bases     = int(cfg.primer3Bases)
    platform = cfg.platform

    seqtype = "iontorrent"
    assert platform.lower() != "illumina", "This module only supports IonTorrent Reads !"
    # set output file prefix
    filePrefixOut = readSet + ".prep"    

    numReads1, numBatches = splitReadFile(readFile1,filePrefixOut,"R1",numCores,deleteLocalFiles)
    # debug check
    if numReads1 == 0:
        raise UserWarning("prep: input read files are empty!")
                                                                
    # run cd-hit to cluster close primer sequences; creates the file {primerFile}.clusters
    primer_trim_ion.cluster_primer_seqs(primerFile)
    
    # set up trimming work to be run in parallel sub-processes, using another python script
    trimScriptIon = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'misc/prep_trim_ion.py')
    workIn = []
    for batchNum in range(numBatches):
        filePrefixBatch = "{}.{:04d}".format(filePrefixOut,batchNum)
        cmd = "python {0} {1} {2} {3} {4} {5} {6} {7} {8} > {6}.log 2>&1 ".format(
            trimScriptIon,cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,
            primer3Bases,filePrefixBatch,primerFile,seqtype)
        workIn.append(cmd)
        
    # run cutadapt and UMI extraction in parallel sub-processes
    print("prep: starting parallel trimming batches")
    pool = ThreadPool(min(numCores,len(workIn)))
    workOut = pool.map(worker, workIn)
    pool.close()
    pool.join()
    print("prep: completed parallel trimming batches")
    
    # make sure all batches of work completed successfully
    for batchNum in range(len(workOut)):
        if not workOut[batchNum]:
            raise Exception("read trimming failed for batch: {:04d}".format(batchNum))
   
    # concatenate the read files back into one file
    for readEnd in ("R1","R2"):
    
        # delete output file if it exists
        readFileOut = "{}.{}.fastq".format(filePrefixOut,readEnd)
        if os.path.isfile(readFileOut):
            os.remove(readFileOut)
      
        # concatenate read file and delete (Linux cat probaby faster than Python line reads)
        for batchNum in range(numBatches):
            readFileIn = "{}.{:04d}.{}.fastq".format(filePrefixOut, batchNum, readEnd)
            cmd = "cat {} >> {} ".format(readFileIn,readFileOut)
            subprocess.check_call(cmd,shell=True)
            os.remove(readFileIn)
   
    # concatenate the log files - note dangerous wildcards here!
    logFileOut = open(filePrefixOut + ".log","w")
    for logFileIn in glob.glob(filePrefixOut + ".0*.log"):
        IN = open(logFileIn,"r")
        logFileOut.write(IN.read())
        IN.close()
        os.remove(logFileIn)
    logFileOut.close()

    # For ion-torrent reads concatenate the umi and primer tag files        
    if seqtype == 'iontorrent':
        OUT1 = open(readSet + ".umi.tag.txt","w")
        OUT2 = open(readSet + ".primer.tag.txt","w")
        OUT3 = open(readSet + ".cutadapt.5.R1.txt","w")
        OUT4 = open(readSet + ".cutadapt.3.R1.txt","w")
        
        for umiTagFileTemp in sorted(glob.glob(filePrefixOut + "*.umi.tag.txt")):
            IN1 =  open(umiTagFileTemp,"r")
            OUT1.write(IN1.read()) # reading into memory
            IN1.close()
            os.remove(umiTagFileTemp)
        for primerTagFileTemp in sorted(glob.glob(filePrefixOut + "*.primer.tag.txt")):
            IN2 = open(primerTagFileTemp,"r")
            OUT2.write(IN2.read()) # reading into memory
            IN2.close()
            os.remove(primerTagFileTemp)
        for cutadapt5FileTemp in sorted(glob.glob(filePrefixOut + "*.cutadapt.5.R1.txt")):
            IN3 = open(cutadapt5FileTemp,"r")
            OUT3.write(IN3.read()) # reading into memory
            IN3.close()
            os.remove(cutadapt5FileTemp)
        for cutadapt3FileTemp in sorted(glob.glob(filePrefixOut + "*.cutadapt.3.R1.txt")):
            IN4 = open(cutadapt3FileTemp,"r")
            OUT4.write(IN4.read()) # reading into memory
            IN4.close()
            os.remove(cutadapt3FileTemp)
            
        OUT1.close()
        OUT2.close()
        OUT3.close()
        OUT4.close()
        
    # aggregate summary read count files - for some trim scripts these files contain important read count metrics
    output = []
    firstFile = True
    for sumFileIn in glob.glob(filePrefixOut + ".0*.summary.txt"):
        iLine = 0
        for line in open(sumFileIn,"r"):
            metricVal,metricName = line.strip().split("\t")
            metricVal = int(metricVal)
            if firstFile:
                output.append([metricVal,metricName])
            else:
                output[iLine][0] += metricVal
            iLine += 1
        firstFile = False
        os.remove(sumFileIn)
    sumFileOut = open(filePrefixOut + ".summary.txt","w")
    for row in output:
        sumFileOut.write("\t".join((str(x) for x in row)))
        sumFileOut.write("\n")
    sumFileOut.close()
