import os
import subprocess
import copy

import pysam

#-------------------------------------------------------------------------------------------------------
# Ion single-end reads - align single-end FASTQ using TMAP
#-------------------------------------------------------------------------------------------------------
def alignToGenomeIon(cfg,readFileIn1,bamFileOut):
    # get some parameters from config
    readSet = cfg.readSet
    numCpus = cfg.numCores
    samtoolsDir = cfg.samtoolsDir
    samtoolsMem = cfg.samtoolsMem
    torrentBinDir = cfg.torrentBinDir
    torrentGenomeFile = cfg.genomeFile
    tmap = os.path.join(torrentBinDir,"tmap")
    tagNameUmiSeq    = cfg.tagNameUmiSeq
    tagNamePrimer    = cfg.tagNamePrimer
    tagNamePrimerErr = cfg.tagNamePrimerErr
    
    # align full-length reads to reference genome using TMAP
    cmd = "{} mapall -n {} -r {} -f {} -v -Y -u --prefix-exclude 5 -o 2 stage1 map4 > ".format(tmap,numCpus,readFileIn1,torrentGenomeFile) \
    + readSet + ".temp.bam 2> " \
    + readSet + ".align.tmap.log "
    subprocess.check_call(cmd, shell=True)

    # add tag to bam
    bamIn = readSet + ".temp.bam"
    bamOut = readSet + ".temp1.bam"
    addBamTags(bamIn,bamOut,readSet,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr)

    # add a fake reverse compliment read alignment (i.e. simulate paired-end primer-side read) for use in downstream code
    bamIn  = pysam.Samfile(readSet + ".temp1.bam", "rb")
    bamOut = pysam.AlignmentFile(bamFileOut, "wb", template=bamIn)
    for read1 in bamIn:
        read1.is_paired = True
        read1.is_read1 = False
        read1.is_read2 = True
        read2 = copy.deepcopy(read1)
        read2.is_read2 = False
        read2.is_read1 = True

        if not read1.is_unmapped:
            read2.is_reverse = not read1.is_reverse
            read1.mate_is_reverse = read2.is_reverse
            read2.mate_is_reverse = read1.is_reverse
            read1.mate_is_unmapped = False
            read2.mate_is_unmapped = False
        bamOut.write(read2)
        bamOut.write(read1)
    bamIn.close()
    bamOut.close()

    # delete unneeded bam files
    os.remove(readSet + ".temp.bam")
    os.remove(readSet + ".temp1.bam")

    # sort by locus for IGV viewing, and for mtMerge with hashing by chromosome
    cmd = samtoolsDir + "samtools sort -m " + samtoolsMem + " -@" + numCpus \
    + " -T " + readSet \
    + " -o " + readSet + ".align.sorted.bam " \
             + readSet + ".align.bam " \
    + "> "   + readSet + ".align.sort.log 2>&1 "
    subprocess.check_call(cmd, shell=True)

    # make BAM index for IGV
    cmd = samtoolsDir + "samtools index " + readSet + ".align.sorted.bam "
    subprocess.check_call(cmd, shell=True)

def addBamTags(bamIn,bamOut,readSet,tagNameUmiSeq,tagNamePrimer,tagNamePrimerErr):
    ''' Add tags to bam
    :param str bamIn : The input bam file to read
    :param str bamOut : The output bam file to write
    :param str readSet : The sample name
    '''
    # create a dict of read id -> tag
    umi_dict = {}
    with open(readSet +  ".umi.tag.txt","r") as IN:
        for line in IN:
            read_id, umi, umi_qual = line.strip('\n').split('\t')
            umi_dict[read_id] = [umi,None,None]
    with open(readSet + ".primer.tag.txt","r") as IN:
        for line in IN:
            read_id, primer, primer_err = line.strip('\n').split('\t')
            umi_dict[read_id][1] = primer
            umi_dict[read_id][2] = primer_err
    print "\nDone creating readID -> (umi,primer,primer_err)  dict\n"            

    with pysam.AlignmentFile(bamIn,"rb") as IN, pysam.AlignmentFile(bamOut,"wb", template=IN) as OUT:
        for read1 in IN:
            temp_tags = read1.tags
            umi_tag = umi_dict["@"+read1.qname][0]
            primer_tag = umi_dict["@"+read1.qname][1]
            primer_err_tag = umi_dict["@"+read1.qname][2]
            temp_tags.append((tagNameUmiSeq,umi_tag))
            temp_tags.append((tagNamePrimer,primer_tag))
            temp_tags.append((tagNamePrimerErr,primer_err_tag))
            read1.tags = tuple(temp_tags)
            OUT.write(read1)
