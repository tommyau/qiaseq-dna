import os
# our modules
from prep import runShellCommand

#-------------------------------------------------------------------------------------
def trim_illumina(filePrefix,cutadaptDir,tagNameUmiSeq):
    ''' Trim and tag synthetic regions from QIASeq DNA Illumina Paired End Reads
    :param str filePrefix: The readSet name
    :param str cutadaptDir: Path to cutadapt executable directory
    :param str tagNameUmiSeq: Tag name for UMI sequence in the SAM/BAM file
    :returns final trimmed R1 and R2 fastq file paths
    :rtype tuple (R1_fastq,R2_fastq)
    '''
    # split read set name from input file directory (NOTE: ".R1.fastq" and ".R2.fastq" are required file suffixes here!)
    dirIn, filePrefix = os.path.split(filePrefix)
    if len(dirIn) > 0:
        dirIn = dirIn + "/"
     
    # trim R1 reads 3' end (gets R2 12 bp barcode, 11-mer common, and ILMN adapter not cut by ILMN software) (primer side)
    cmd  = cutadaptDir + "cutadapt -e 0.18 -O 3 " \
         + "-a AGGACTCCAAT -n 1 " \
         + "-o " + filePrefix + ".temp0.R1.fastq " \
         + dirIn + filePrefix +       ".R1.fastq "
    runShellCommand(cmd)
    
    # trim R2 reads 3' end (gets R1 custom sequencing adapter region not cut by ILMN software - customer sequencing primer settings usually wrong) (barcode side)
    cmd  = cutadaptDir + "cutadapt -e 0.18 -O 3 " \
         + "-a CAAAACGCAATACTGTACATT -n 1 " \
         + "-o " + filePrefix + ".temp0.R2.fastq " \
         + dirIn + filePrefix +       ".R2.fastq "
    runShellCommand(cmd)
       
    # open output
    fileOut1 = open(filePrefix + ".temp1.R1.fastq", "w")
    fileOut2 = open(filePrefix + ".temp1.R2.fastq", "w")
    fileIn1  = open(filePrefix + ".temp0.R1.fastq", "r")
    fileIn2  = open(filePrefix + ".temp0.R2.fastq", "r")
 
    # extract barcode at beginning of R1 - just take first 12 bp
    numReadPairsTotal = 0
    numReadPairsDropped = 0
    linesR1 = []
    linesR2 = []
    lineIdx = 0
    for line in fileIn1:
        linesR1.append(line.strip())
        line = fileIn2.readline()
        linesR2.append(line.strip())
        lineIdx += 1
  
        # read all four lines of R1 and R2 into memory
        if lineIdx == 4:
            lineIdx = 0
            numReadPairsTotal += 1
         
            # skip reads too short
            if len(linesR1[1]) < 40 or len(linesR2[1]) < 40:
                numReadPairsDropped += 1
                del(linesR1[:])
                del(linesR2[:])
                continue
             
            # get barcode, trim R2 - 12 bp barcode, 11 bp univeral
            line = linesR2[1]
            if not line.startswith("N"):
                umiSeq = line[0:12]
                linesR2[1] = linesR2[1][23:]
                linesR2[3] = linesR2[3][23:]
            else:
                umiSeq = line[1:13]
                linesR2[1] = linesR2[1][24:]
                linesR2[3] = linesR2[3][24:]
    
            # make a UMI tag
            umiTag = tagNameUmiSeq + ":Z:" + umiSeq
               
            # add barcode to R1 id line
            line = linesR1[0]
            idx = line.find(" ")
            readId = line[1:idx]
            linesR1[0] = "@" + readId + " " + umiTag
            
            # debug check on R2 being in sync with R1
            line = linesR2[0]
            idx = line.find(" ")
            readId_ = line[1:idx]
            if readId_ != readId:
                raise Exception("R1/R2 no synchronized")
             
            # add barcode to R2 id line
            linesR2[0] = linesR1[0]
            
            # write output
            for i in range(4):
                fileOut1.write(linesR1[i] + "\n")
                fileOut2.write(linesR2[i] + "\n")
    
            # clear for next read           
            del(linesR1[:])
            del(linesR2[:])
      
    fileOut1.close()
    fileOut2.close()
    
    # write summary file 
    fileOut = open(filePrefix + ".prep_trim.summary.txt", "w")
    fileOut.write("{}\tread fragments total\n".format(numReadPairsTotal))
    fileOut.write("{}\tread fragments dropped, less than 40 bp\n".format(numReadPairsDropped))
    fileOut.close()

    return (filePrefix + ".temp1.R1.fastq",filePrefix+".temp1.R2.fastq")

#-------------------------------------------------------------------------------------------------------
# trimming using recognition of universal sequence on both sides, and extract duplex tag TT or CC
#-------------------------------------------------------------------------------------------------------
def trim_illumina_duplex(filePrefix,cutadaptDir,tagNameUmiSeq):
    ''' Trim and tag synthetic regions from QIASeq Dna Duplex Paired End Reads
    :param str filePrefix: The readSet name
    :param str cutadaptDir: Path to cutadapt executable directory
    :param str tagNameUmiSeq: Tag name for UMI sequence in the SAM/BAM file
    :returns final trimmed R1 and R2 fastq file paths
    :rtype tuple (R1_fastq,R2_fastq)
    '''
    # Note : R1 and R2 are swapped for the trimming section (have an option here to signal when to do this)

    # split read set name from input file directory (NOTE: ".R1.fastq" and ".R2.fastq" are required file suffixes here!)
    dirIn, filePrefix = os.path.split(filePrefix)
    if len(dirIn) > 0:
        dirIn = dirIn + "/"

    # trim R1 reads (primer side), 3' end (also gets 12 bp barcode and ILMN adapter not cut by ILMN software)
    cmd  = cutadaptDir + "cutadapt -e 0.18 -O 3 " \
         + "-a AGGACTCCAAT -n 1 " \
         + "-o " + filePrefix + ".temp0.R2.fastq " \
         + dirIn + filePrefix + ".R1.fastq "
    runShellCommand(cmd)

    # trim R2 reads (barcode side), 3' end (gets ILMN adapter not cut by ILMN software - customer sequencing primer settings usually wrong)
    cmd  = cutadaptDir + "cutadapt -e 0.18 -O 3 " \
         +  "-a CAAAACGCAATACTGTACATT -n 1 " \
         + "-o " + filePrefix + ".temp0.R1.fastq " \
         + dirIn + filePrefix + ".R2.fastq "         
    runShellCommand(cmd)

    # trim 5' end of R2 , discard read pairs where the universal is not present on the 5' end of R2 
    cmd = cutadaptDir + "cutadapt -e 0.18 -O 18 --discard-untrimmed --minimum-length 35 " \
        + "-u 12 -g ^TTCTGAGCGAYYATAGGAGTCCT " \
        + "-o {0}.temp1.R1.fastq -p {0}.temp1.R2.fastq {0}.temp0.R1.fastq {0}.temp0.R2.fastq ".format(filePrefix)
    runShellCommand(cmd)

    # init counters
    numReadPairsTotal = 0
    numReadPairsDropped = 0
    numReadPairsTT = 0
    numReadPairsCC = 0
    numReadPairsNN = 0

    # open input and output files
    fileout1 = open(filePrefix + ".temp2.R2.fastq","w")
    fileout2 = open(filePrefix + ".temp2.R1.fastq","w")
    filein1  = open(filePrefix + ".temp1.R1.fastq","r")
    filein2  = open(filePrefix + ".temp1.R2.fastq","r")
    filein3  = open(filePrefixOut + ".temp0.R1.fastq","r")  # needed to pull barcode

    # loop over all R1 reads
    lineIdx = 0
    for line1 in filein1:
        line2 = filein2.readline()

        # echo lines 2,3,4
        if lineIdx > 0:
            fileout1.write(line1)
            fileout2.write(line2)

        # add read id and duplex tag to line 1 readId
        else:
            idx = line1.find(" ")
            readId1 = line1[:idx]
            readId2 = line2[:idx]
            assert readId2 == readId1,"R1 and R2 do not have the same name"

            # get the read sequence from temp0.R2.fastq file, which contains all reads, all 5' bases
            line3 = None
            while True:
                line3 = filein3.readline()
                numReadPairsTotal += 1
                i = line3.find(" ")
                readId3 = line3[:i]
                line3 = filein3.readline()
                filein3.readline() # discard line 3
                filein3.readline() # discard line 4
                if readId3 == readId1:
                    break
                else:
                    numReadPairsDropped += 1
            if line3 == None:
                raise Exception()

            # extract barcode and duplex tag (rough initial cut)
            umi = line3[0:12]
            duplexTagRegion = line3[21:25]
            if duplexTagRegion.find("CC") >= 0:
                duplexTag = "CC"
                numReadPairsCC += 1
            elif duplexTagRegion.find("TT") >= 0:
                duplexTag = "TT"
                numReadPairsTT += 1
            else:
                duplexTag = "NN"
                numReadPairsNN += 1

            # add the duplex tag to read id and the umi tag  after the space so that bwa can create a seperate tag for it
            # note : bwa mem -C does not seem to create 2 tags
            umiTag = tagNameUmiSeq + ":Z:" + umi
            duplexTag = "DU:Z:" + duplexTag

            line1 = readId1 + "\t" +  duplexTag + "\t" + umiTag
            line2 = readId1 + "\t" +  duplexTag + "\t" + umiTag
            fileout1.write(line1+"\n")
            fileout2.write(line2+"\n")

        # prepare next iteration
        lineIdx += 1
        if lineIdx == 4:
            lineIdx = 0

    # close files
    fileout1.close()
    fileout2.close()
    filein1.close()
    filein2.close()
    filein3.close()

    # write summary file
    fileout = open(filePrefix + ".prep_trim.summary.txt", "w")
    fileout.write("{}\tread fragments total\n".format(numReadPairsTotal))
    fileout.write("{}\tread fragments dropped, less than 40 bp\n".format(numReadPairsDropped))
    fileout.write("{}\tread fragments with duplex tag CC\n".format(numReadPairsCC))
    fileout.write("{}\tread fragments with duplex tag TT\n".format(numReadPairsTT))
    fileout.write("{}\tread fragments with duplex tag NN\n".format(numReadPairsNN))
    fileout.close()

    return (filePrefix + ".temp2.R1.fastq",filePrefix + ".temp2.R2.fastq")
