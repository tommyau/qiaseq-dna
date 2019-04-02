import sys

import numpy as np


def run(cfg):
    '''
    '''
    print('starting fraglen_by_rpu...')
    
    # input
    readset                 = cfg.readSet
    umi_mark_for_sum_primer = readset + '.umi_mark.for.sum.primer.txt'

    # output
    outfile = readset + '.fraglen_by_rpu.distrib.txt'
    out_header = "|".join(
        ("read set", "read frags per umi", "num umis", "num read frags" , "mean frag len", \
         "25th percentile frag len", "50th percentile frag len", "75th percentile frag len"))
    
    
    umis_with_1_readfrag = []
    umis_with_2_readfrag = []
    umis_with_3_readfrag = []
    umis_1rpu = 0
    read_frags_1rpu = 0
    umis_2rpu = 0
    read_frags_2rpu = 0
    umis_3rpu = 0
    read_frags_3rpu = 0

    with open(umi_mark_for_sum_primer, "r") as IN:
        for line in IN:
            chrom,strand,umi_loc,umi,num_reads,num_alignments,mt_read_idx,is_resample,frag_len,primer,primer_loc5 = line.strip("\n").split("|")

            loc3 = int(primer_loc5) + len(primer) - 1 if strand == "0" else int(primer_loc5) - len(primer)
            molecule = "-".join([chrom, strand, umi_loc, umi])
            num_reads = int(num_reads)

            if num_reads == 1: # this is read frags
                umis_with_1_readfrag.append(int(frag_len))
                umis_1rpu += 1
                read_frags_1rpu += num_reads


            elif num_reads == 2: # this is read frags
                umis_with_2_readfrag.append(int(frag_len))
                umis_2rpu += 1
                read_frags_2rpu += num_reads

            elif num_reads == 3: # this is read frags
                umis_with_3_readfrag.append(int(frag_len))
                umis_3rpu += 1
                read_frags_3rpu += num_reads

    umis_with_1_readfrag = np.array(umis_with_1_readfrag)
    umis_with_2_readfrag = np.array(umis_with_2_readfrag)
    umis_with_3_readfrag = np.array(umis_with_3_readfrag)

    num_umis = [umis_1rpu, umis_2rpu, umis_3rpu]
    num_reads = [read_frags_1rpu, read_frags_2rpu, read_frags_3rpu]

    with open(outfile, "w") as OUT:
        OUT.write(out_header)
        OUT.write("\n")
        for i,arr in enumerate([umis_with_1_readfrag,umis_with_2_readfrag,umis_with_3_readfrag]):
            percentile_25 = np.percentile(arr,25)
            percentile_50 = np.percentile(arr,50)
            percentile_75 = np.percentile(arr,75)
            mean = round(np.mean(arr),2)

            umis = num_umis[i]
            reads = num_reads[i]

            temp = [readset,i+1,umis,reads,mean,percentile_25,percentile_50,percentile_75]
            OUT.write("|".join([str(x) for x in temp]))
            OUT.write("\n")


if __name__ == '__main__':
    readset   = sys.argv[1]
    paramfile = sys.argv[2]
    import core.run_config
    cfg = core.run_config.run(readset, paramfile)
    run(cfg)
        
