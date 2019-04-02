import os
import subprocess
import sys
import collections
import itertools
import edlib
import cutadapt.adapters
import time
import cPickle as pickle

# our modules
from core.umi_filter import reverseComplement

'''
To do :
1.) Add Proper Command Line Parser
2.) Output some metrics like num_bases trimmed, num_reads trimmed etc.
3.) Add parameterization for k-mer index, min error rate, etc.
'''

def cluster_primer_seqs(primer_file):
    ''' Cluster similar primer sequences for use later in primer trimming
    :param str primerFile: Path to the primerfile for this readset
    '''
    # create a fasta of the primers
    cmd1 = (
       """ i=0; while read chrom pos strand primer; do echo ">"$i; echo $primer|tr -d '\r'; """
       """ i=$(($i+1)); done < {primerfile} > {primerfile}.fasta; """.format(primerfile=primer_file)
       )
    print cmd1
    subprocess.check_call(cmd1,shell=True)
    # run cd-hit
    cmd2 = "/srv/qgen/bin/downloads/cd-hit-v4.6.8-2017-1208/cd-hit-est -i {primerfile}.fasta -o {primerfile}.clusters.temp".format(
       primerfile=primer_file)
    subprocess.check_call(cmd2,shell=True)
    # parse output to usable format
    parse_cdhit(primer_file,primer_file+'.clusters.temp.clstr',primer_file+'.clusters')


def parse_cdhit(primer_file,cdhit_out,simple_cdhit_out):
    ''' Parse cd-hit output to usabel format
    :param str primer_file: Path to the the primerfile for this readset
    :param str cdhit_out: Path to output file from cd-hit
    :param str simple_cdhit_out: Path to write the parsed results to
    '''
    primers = []
    i=0
    with open(primer_file,'r') as IN:
        for line in IN:
            primers.append(line.strip('\n\r').split('\t')[-1])

    cluster_info = collections.defaultdict(list)
    with open(cdhit_out,'r') as IN:
        for line in IN:
            if line.startswith('>'):
                cluster_num = int(line.strip('\n').split(' ')[1])
                new_cluster = True
            else:
                temp = line.strip('\n').split('\t')[1].split(',')[1].strip().split('...')[0].strip('>')
                primer = primers[int(temp)]
                cluster_info[str(cluster_num)].append(primer)

    with open(simple_cdhit_out,'w') as OUT:
        for cluster in cluster_info:
            if len(cluster_info[cluster]) > 1:
                OUT.write(','.join(cluster_info[cluster])+'\n')


def create_primer_search_datastruct(primer_file,primer_file_clusters,cache=False,cache_file=None):
    ''' Create datastructures used for searching primers
    Overlapping 8-mer index for each primer
    Information for each primer including closely clustered primers from cd-hit

    :param str primer_file: The PrimerFile for this readset
    :param str primer_file_clusters : The Clustering output from cd-hit
    :param bool cache: Whether to cache the results to a file
    :param str cache_file: If cache is True , file path to cache to
    :rtype If cache is False, a tuple of : (defaultdict(set),defaultdict(list),dict)
    '''
    primer_kmer = collections.defaultdict(set)
    primers = collections.defaultdict(list)
    primers_cutadapt = {}

    with open(primer_file,'r') as IN:
        for line in IN:
            contents = line.strip('\n\r').split('\t')
            primer = contents[-1]
            if primer in primers:
                raise Exception("Duplicate Primer Encountered : %s"%line)
            contents.append(len(primer))
            contents.append([]) # place-holder for clustered primer sequences
            primers[primer] = contents

            # create k-mer index
            k=8
            kmers = set(''.join(itertools.islice(primer,i,i+k)) for i in range(len(primer)+1-k))
            for oligo in kmers:
                primer_kmer[oligo].add(primer)
            # create object for cutadapt
            revcomp = reverseComplement(primer)
            primers_cutadapt[primer] = cutadapt.adapters.Adapter(sequence=revcomp,where=cutadapt.adapters.BACK,max_error_rate=0.1,min_overlap=3)

    # add close primer sequences for each primer based on cdhit results
    with open(primer_file_clusters,'r') as IN:
        for line in IN:
            temp = line.strip('\n').split(',')
            for p1,p2 in itertools.product(temp, repeat=2):
                assert p1 in primers and p2 in primers,"Primer(s) from cd-hit not in PrimerFile !!!"
                primers[p1][-1].append(p2)
                primers[p2][-1].append(p1)
    if cache:
        with open(cache_file,"wb") as OUT:
            pickle.dump((primer_kmer,primers,primers_cutadapt),OUT)
    else:
        return (primer_kmer,primers,primers_cutadapt)

def trim_primer(primer_datastruct,R1):
    ''' Trim appropriate Primer from the R1 read
    :param tuple primer_datastruct: Object returned by the function create_primer_search_datastruct
    :param R1 : The R1 read sequence
    :returns : The best primer hit and the trimming site
    :rtype tuple
    '''
    k = 8
    best_editdist = None
    best_score = None
    best_primer = None
    best_plen = None
    primer_kmer , primers, primers_cutadapt = primer_datastruct

    candidates = set()
    for oligo in set(''.join(itertools.islice(R1[0:30],i,i+8)) for i in range(len(R1[0:30])+1-k)): # all 8-mers of the first 30 bp of the read(length will be truncated to read length if read is less than 30 bp)
        if oligo in primer_kmer: # check if the 8-mer is in the index
            for c in primer_kmer[oligo]: # iterate over all primers corresponding to this 8-mer
                candidates.add(c)
                for similar_primer in primers[c][-1]: # iterate over all similar primers to this primer
                    candidates.add(similar_primer)

    if len(candidates) == 0: # no hits in the index, exhaustive search over all primers
        candidates = primers

    num_candidates = len(candidates)
    for p in candidates:

        p_len = primers[p][-2]

        if R1[0:p_len] == p and num_candidates == 0: # exact match and no other primers to check
            return (p_len+1,p,'0')

        else:
            alignment = edlib.align(p,R1[0:p_len+5],mode="SHW")
            editdist = alignment['editDistance']
            score =  float(editdist)/p_len
            if best_score == None or score < best_score:
                temp = alignment
                best_primer = p
                best_score = score
                best_editdist = editdist
                best_plen = p_len
            elif score == best_score:
                if best_plen < p_len: # when same score with same primer length; the first hit will be chosen
                    temp = alignment
                    best_primer = p
                    best_score = score
                    best_editdist = editdist
                    best_plen = p_len

    assert best_score != None, "Bug in logic ! Primer could not be scored correctly !"

    if best_score <= 0.10:
        return (temp['locations'][-1][1]+1,best_primer,str(best_editdist))
    else:
        return (-1,"-1","-1")

def iterate_fastq(R1_fastq,R2_fastq):
    '''
    '''
    with open(R1_fastq,'r') as IN1 , open(R2_fastq,'r') as IN2:
        while True:
            R1_info =  (IN1.next().strip('\n').strip('\n'),IN1.next().strip('\n'),IN1.next().strip('\n'),IN1.next().strip('\n'))
            R2_info =  (IN2.next().strip('\n'),IN2.next().strip('\n'),IN2.next().strip('\n'),IN2.next().strip('\n'))
            yield (R1_info,R2_info)



def main(R1_fastq,R2_fastq,R1_fastq_out,R2_fastq_out,primer_file,primer_file_clusters,primer_3_bases,primer_tag_name,primer_err_tag_name,update_read_id=True,out_tag_file=None,load_cache=False,cache_file=None):
    ''' Main function
    :param str R1_fastq: Path to Input R1 fastq file
    :param str R2_fastq: Path to Input R2 fastq file
    :param str R1_fastq_out: Path to Output R1 fastq file
    :param str R2_fastq_out: Path to Output R2 fastq file
    :param str primer_file: Path to primer file for this readset
    :param str primer_file_clusters: Path to the cdhit clustering output of the primers
    :param int primer_3_bases: Number of bases to keep on the 3' end of the primers
    :param str primer_tag_name: Tag name for storing the primer in the bam/sam file
    :param str primer_err_tag_name: Tag Name for storing the edit distance of the primer match in the bam/sam file
    :param bool update_read_id: Whether to add primer info tags to read id. This is False for IonTorrent reads since tmap cannot tag bam files with these
    :param str out_tag_file : Output file for primer info tags if update_read_id is False. None by default  
    :param bool load_cache: Whether to load the primer search datastructure (useful when splitting a fastq and running in parallel)
    :param str cache_file: If load_cache is True, file path to load from  
    '''
    # counters
    num_R1=0
    trimmed_R1=0
    trimmed_R2= 0

    # primer search datastruct
    if load_cache: # pickle and dill do not seem to be working. will debug this later. do not use cache for now.
        raise UserWarning("Feature not implemented yet, cannot serialize primer data structure objects.")
        with open(cache_file,"rb") as IN:
            primer_datastruct = pickle.load(IN)
    else:
        primer_datastruct = create_primer_search_datastruct(primer_file,primer_file_clusters,cache=False)
    primer_kmer,primers,primers_cutadapt = primer_datastruct

    def reformat_readid(read_id,primer_info,primer_err,update_read_id=True):
        ''' Reformat read id with primer tags
        '''
        if update_read_id:
            primer_info_tag = primer_tag_name + ":Z:" + primer_info
            primer_err_tag = primer_err_tag_name + ":Z:" + primer_err
        else:
            primer_info_tag = primer_info
            primer_err_tag = primer_err
        return (read_id + "\t" + primer_info_tag + "\t" + primer_err_tag)

    if not update_read_id:
        OUT3 = open(out_tag_file,"w")       
        
    with open(R1_fastq_out,'w') as OUT1, open(R2_fastq_out,'w') as OUT2:
        for R1_info, R2_info in iterate_fastq(R1_fastq,R2_fastq):
            R1_id,R1_seq,R1_t,R1_qual = R1_info
            R2_id,R2_seq,R2_t,R2_qual = R2_info
            temp_id = R1_id.split(" ")[0]
            trim_pos,primer,primer_err = trim_primer(primer_datastruct,R1_seq)

            if trim_pos == -1: # No primer found
                # re-format read id
                if update_read_id:
                    R1_id = reformat_readid(R1_id,primer,primer_err)
                    R2_id = reformat_readid(R2_id,primer,primer_err)
                else:
                    OUT3.write(reformat_readid(R1_id,primer,primer_err,update_read_id)+"\n") # for ion torrent reads                    
                # write to output fastq
                OUT1.write(R1_id + "\n" + R1_seq + "\n" + R1_t + "\n" + R1_qual + "\n")
                OUT2.write(R2_id + "\n" + R2_seq + "\n" + R2_t + "\n" + R2_qual + "\n")

            else: # found primer
                chrom,pos,strand,seq,p_len,similar_primers = primers[primer]
                primer_info = chrom+"-"+strand+"-"+pos+"-"+str(len(primer))
                trimmed_R1+=1
                # re-format read id
                if update_read_id:
                    R1_id = reformat_readid(R1_id,primer_info,primer_err)
                    R2_id = reformat_readid(R2_id,primer_info,primer_err)
                else:
                    OUT3.write(reformat_readid(R1_id,primer_info,primer_err,update_read_id)+"\n") # for ion torrent reads
                # trim R1
                if primer_3_bases ==  -1 or primer_3_bases > len(R1_seq): # keep R1 and R1_qual to be as is
                    pass
                elif primer_3_bases == 0: # trim all bases belonging to the primer
                    R1_seq = R1_seq[trim_pos:]
                    R1_qual = R1_qual[trim_pos:]
                else: # keep said bases belonging to the primer
                    R1_seq = R1_seq[trim_pos - primer_3_bases:]
                    R1_qual = R1_qual[trim_pos - primer_3_bases:]

                OUT1.write(R1_id + "\n" + R1_seq + "\n" + R1_t + "\n" + R1_qual + "\n")

                # trim R2
                read_seq = cutadapt.seqio.Sequence(name=R2_id,sequence=R2_seq)
                match = primers_cutadapt[primer].match_to(read=read_seq)
                if match == None: # could not match primer
                    OUT2.write(R2_id + "\n" + R2_seq + "\n" + R2_t + "\n" + R2_qual + "\n")
                else:
                    trimmed_R2+=1
                    if primer_3_bases == -1 or primer_3_bases > len(R2_seq): # keep R1 and R1_qual to be as is
                        pass
                    elif primer_3_bases == 0: # trim all bases belonging to the primer
                        R2_seq = R2_seq[0:match.rstart]
                        R2_qual = R2_qual[0:match.rstart]
                    else: # trim while keeping said bases from the primer; note R2 will be truncated to its length if primer region found is less than primer_3_bases
                        R2_seq = R2_seq[0:match.rstart+primer_3_bases]
                        R2_qual = R2_qual[0:match.rstart+primer_3_bases]

                    OUT2.write(R2_id+ "\n" + R2_seq+"\n" + R2_t + "\n" + R2_qual + "\n")

            num_R1+=1

    print "Total R1 Reads: {}".format(num_R1)
    print "R1 Reads Trimmed: {}".format(trimmed_R1)
    print "Total R2 Reads: {}".format(num_R1)
    print "R2 Reads Trimmed: {}".format(trimmed_R2)
    if not update_read_id:
        OUT3.close()

if __name__ == '__main__':
    assert len(sys.argv) == 10, "Incorrect command line params specified !"
    main(*sys.argv[1:])
