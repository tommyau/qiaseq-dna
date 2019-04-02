import sys
import string
import subprocess
import os

# our modules
import primer_trim_ion

#-------------------------------------------------------------------------------------
def trimIon(filePrefix,cutadaptDir,tagNameUmiSeq):
    ''' Trim and tag synthetic regions from QIASeq Dna IonTorrent Reads
    Please note that this function is tied to the core/prep.py and core/prep_trim.py modules
    The input and output files are split and merged in those modules.

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

    # look for universal seq at 3' end of read, drop read if not found (these are long fragments)
    cmd = cutadaptDir + "cutadapt -e 0.18 -O 17 --discard-untrimmed " \
        + "-a CAAAACGCAATACTGTACATT -n 1 " \
        + "--info-file " + filePrefix + ".cutadapt.3.R1.txt " \
        + "-o " + filePrefix + ".temp0.R1.fastq " \
        + dirIn + filePrefix + ".R1.fastq" \
        + " > " + filePrefix + ".cutadapt.3.R1.log 2>&1 "
    subprocess.check_call(cmd, shell=True)

    # read original read count from cutadapt log
    numReadsTotal = None
    for line in open(filePrefix + ".cutadapt.3.R1.log","r"):
        if line.startswith("Total reads processed:"):
            (key,val) = line.strip().split(":")
            numReadsTotal = int(val.strip().replace(",",""))
            break
    os.remove(filePrefix + ".cutadapt.3.R1.log")

    # pull first 12 bp off, move to header
    fileout = open(filePrefix + ".temp1.R1.fastq","w")
    fileoutTag = open(filePrefix + ".umi.tag.txt","w")
    numReadsWithAdapter3 = 0
    numReadsDroppedTooShort = 0
    lines = []
    lineIdx = 0
    for line in open(filePrefix + ".temp0.R1.fastq","r"):
        lines.append(line.strip())
        lineIdx += 1

        # all four lines of R1 in memory
        if lineIdx == 4:
            lineIdx = 0
            numReadsWithAdapter3 += 1

            # skip reads too short
            if len(lines[1]) < 40:
                numReadsDroppedTooShort += 1
                del(lines[:])
                continue

            # get barcode, trim R1 - 12 bp barcode
            barcode  = lines[1][0:12]
            barcodeQ = lines[3][0:12]
            lines[1] = lines[1][12:]
            lines[3] = lines[3][12:]

            # not adding umi to read id as tmap cannot add bam tags from read id comment
            line = lines[0]

            # write output fastq
            for i in range(4):
                fileout.write(lines[i] + "\n")

            # write umi to file, add as tag to bam later
            fileoutTag.write(lines[0]+"\t"+barcode+"\t"+barcodeQ+"\n")

            # clear for next read
            del(lines[:])
    fileout.close()
    numReadsDroppedNoAdapter3 = numReadsTotal - numReadsWithAdapter3

    # trim 11-mer universal adapter from 5' end
    cmd = cutadaptDir + "cutadapt -e 0.18 -O 9 --discard-untrimmed " \
        + "-g ^ATTGGAGTCCT -n 1 " \
        + "--info-file " + filePrefix + ".cutadapt.5.R1.txt " \
        + "-o " + filePrefix + ".temp2.R1.fastq " \
        + filePrefix + ".temp1.R1.fastq" \
        + " > " + filePrefix + ".cutadapt.5.R1.log 2>&1 "
    subprocess.check_call(cmd, shell=True)

    # get read count from cutadapt log
    numReadsFinal = 0
    for line in open(filePrefix + ".cutadapt.5.R1.log","r"):
        if line.startswith("Reads with adapters:"):
            (key,val) = line.strip().split(":")
            (val,foo) = val.split("(")
            numReadsFinal = int(val.strip().replace(",",""))
            break
    numReadsDroppedNoAdapter5 = numReadsTotal - numReadsDroppedNoAdapter3 - numReadsDroppedTooShort - numReadsFinal
    os.remove(filePrefix + ".cutadapt.5.R1.log")

    # write summary file
    fileout = open(filePrefix + ".align.summary.txt", "w")
    fileout.write("{}\tread fragments total\n".format(numReadsTotal))
    fileout.write("{}\tread fragments dropped, not full-length\n".format(numReadsDroppedNoAdapter3))
    fileout.write("{}\tread fragments dropped, less than 40 bp\n".format(numReadsDroppedTooShort))
    fileout.write("{}\tread fragments dropped, universal 11-mer not found\n".format(numReadsDroppedNoAdapter5))
    fileout.close()

    # set up reverse comlement
    dnaComplementTranslation = string.maketrans("ATGC", "TACG")

    # make fake R2 (primer side) file
    fileout = open(filePrefix + ".temp2.R2.fastq", "w")
    lineIdx = 0
    for line in open(filePrefix + ".temp2.R1.fastq", "r"):
        line = line.strip()
        if lineIdx == 0 or lineIdx == 2:
            fileout.write(line)
        elif lineIdx == 3:
            fileout.write(line[::-1])
        else:
            seq = line[::-1]
            seq = seq.translate(dnaComplementTranslation)
            fileout.write(seq)
        fileout.write("\n")
        lineIdx += 1
        if lineIdx == 4:
            lineIdx = 0
    fileout.close()

    return (filePrefix + ".temp2.R1.fastq",filePrefix + ".temp2.R2.fastq")

#-------------------------------------------------------------------------------------
def run(cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,primer3Bases,filePrefix,primerFile,seqtype):

    print "Trimming synthetic regions , i.e. Adapters, UMI\n"
    fileIn1,fileIn2 = trimIon(filePrefix,cutadaptDir,tagNameUmiSeq)
    print "Trimming Primers\n"
    fileOut1 = filePrefix + ".trimmed.R1.fastq"
    fileOut2 = filePrefix + ".trimmed.R2.fastq"
    primer_trim_ion.main(
        fileIn2,fileIn1,fileOut2,fileOut1,primerFile,primerFile+'.clusters',
        int(primer3Bases),tagNamePrimer,tagNamePrimerErr,update_read_id=False,
        out_tag_file=filePrefix+".primer.tag.txt",load_cache=False,cache_file=None) # swapped R1 and R2
    print "Finished with trimming\n"    
    # delete unneeded temp files
    for i in range(0,3):
        for read in ('R1','R2'):
            file_to_delete = filePrefix + ".temp%i.%s.fastq"%(i,read)
            if os.path.exists(file_to_delete):
                os.remove(file_to_delete)
    # rename the output files - overwrite the input files!
    os.rename(filePrefix + ".trimmed.R1.fastq", filePrefix + ".R1.fastq")
    os.rename(filePrefix + ".trimmed.R2.fastq", filePrefix + ".R2.fastq")

#-------------------------------------------------------------------------------------
if __name__ == "__main__":
    cutadaptDir = sys.argv[1]
    tagNameUmiSeq = sys.argv[2]
    tagNamePrimer = sys.argv[3]
    tagNamePrimerErr = sys.argv[4]
    primer3Bases = sys.argv[5]
    readFilePrefix = sys.argv[6]
    primerFile = sys.argv[7]
    seqtype = sys.argv[8]
    run(cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,primer3Bases,readFilePrefix,primerFile,seqtype)
