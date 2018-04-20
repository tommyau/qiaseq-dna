import os
import sys
from collections import defaultdict

def cluster_primer_seqs(primer_file):
   ''' Cluster similar primer sequences for use later in primer trimming 
   :param str primerFile: Path to the primerfile for this readset 
   '''
   # create a fasta of the primers
   cmd1 = (
      """ i=0; while read chrom pos strand primer; do echo ">"$i; echo $primer; """
      """ i=$(($i+1)); done < {primerfile} > {primerfile}.fasta; """.format(primerfile=primer_file)
      )
   subprocess.check_call(cmd1)
   # run cd-hit
   cmd = "/srv/qgen/bin/downloads/cd-hit-v4.6.8-2017-1208/cd-hit -i {primerfile}.fasta -o {primerfile}.clusters.temp".format(
      primerfile=primer_file)
   subprocess.check_call(cmd)
   # parse output to usable format
   parse_cdhit(primer_file,primer_file+'.clusters.temp',primer_file+'.clusters')


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
            primers.append(line.strip('\n'))

    cluster_info = defaultdict(list)
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
