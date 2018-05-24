import sys
import subprocess
import os
import os.path
# our modules
from prep_trim_options import trim_illumina, trim_illumina_duplex
from misc import trim_iontorrent
import primer_trim

# Need to choose appropriate trimming function for different sequencing types
_trimming_function_ = {
'illumina'            :  trim_illumina, 
'illumina_duplex'     :  trim_illumina_duplex, 
'iontorrent'          :  trim_iontorrent
}

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
def run(cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,primer3Bases,filePrefix,primerFile,seqtype):

    print "Trimming synthetic regions , i.e. Adapters, UMI\n"
    fileIn1,fileIn2 = _trimming_function_[seqtype](filePrefix,cutadaptDir,tagNameUmiSeq)
    print "Trimming Primers\n"
    fileOut1 = filePrefix + ".trimmed.R1.fastq"
    fileOut2 = filePrefix + ".trimmed.R2.fastq"   
    primer_trim.main(fileIn1,fileIn2,fileOut1,fileOut2,primerFile,primerFile+'.clusters',int(primer3Bases),tagNamePrimer,tagNamePrimerErr,load_cache=True,cache_file=primerFile+".kmer.cache")
    print "Finished with trimming\n"    
    # delete unneeded temp files
    for i in range(0,3):
        file_to_delete = filePrefix + ".temp%i.R1.fastq"%i
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
    primerFasta = sys.argv[5]
    primer3Bases = sys.argv[6]
    readFilePrefix = sys.argv[7]
    primerFile = sys.argv[8]
    seqtype = sys.argv[9]
    run(cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,primerFasta,primer3Bases,readFilePrefix,primerFile,seqtype)
