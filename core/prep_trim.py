import sys
import subprocess
import os
import os.path
# our modules
from prep_trim_options import trim_illumina, trim_illumina_duplex
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from misc.process_ion import trim_iontorrent
import primer_trim

# Need to choose appropriate trimming function for different sequencing types
_trimming_function_ = {
'illumina'            :  trim_illumina, 
'illumina_duplex'     :  trim_illumina_duplex, 
'iontorrent'          :  trim_iontorrent
}

#-------------------------------------------------------------------------------------
def run(cutadaptDir,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr,primer3Bases,filePrefix,primerFile,seqtype):

    print "Trimming synthetic regions , i.e. Adapters, UMI\n"
    fileIn1,fileIn2 = _trimming_function_[seqtype](filePrefix,cutadaptDir,tagNameUmiSeq)
    print "Trimming Primers\n"
    fileOut1 = filePrefix + ".trimmed.R1.fastq"
    fileOut2 = filePrefix + ".trimmed.R2.fastq"   
    primer_trim.main(fileIn1,fileIn2,fileOut1,fileOut2,primerFile,primerFile+'.clusters',int(primer3Bases),tagNamePrimer,tagNamePrimerErr,load_cache=False,cache_file=primerFile+".kmer.cache")
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
