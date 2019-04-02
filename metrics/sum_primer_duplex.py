from collections import defaultdict
import sys

import pysam


def run(cfg):
    '''
    '''
    print('sum_primer_duplex starting...')
    
    # input
    readset       = cfg.readSet
    inbam         = readset + '.umi_merge.bam'
    umi_mark_file = readset + '.umi_mark.for.sum.primer.txt'

    # output
    outfile       = readset + '.sum.primer.duplex.txt'
    # output metric file header
    out_header = "|".join(
        ["read set", "primer", "strand", "chrom", "loc5", "loc3", "Total UMIs", \
         "Total read frags","CC UMIs","CC read frags","TT UMIs", \
         "TT read frags", "NN UMIs", "NN read frags", \
         "NN Only UMIs", "Duplex UMIs (1 read frag CC and TT)"])

    # store duplex tag for each UMI
    IN = pysam.AlignmentFile(inbam,"rb")
    duplex_by_umi = defaultdict(lambda:defaultdict(int))    
    for read in IN:
        umi = read.get_tag(cfg.tagNameUmi)
        duplex_id = read.get_tag(cfg.tagNameDuplex) 
        duplex_by_umi[umi][duplex_id] += 1
        
    IN.close()

    primer_metrics = defaultdict(lambda:defaultdict(int))
    primers = set()
    primer_info = {}

    # store duplex UMIs for each primer
    with open(umi_mark_file, "r") as IN:
        for line in IN:
            chrom, strand, umi_loc, umi, num_read_frags, num_alignments, mt_read_idx, \
                is_resample,frag_len,primer,primer_loc5 = line.strip().split("|")            
            primers.add(primer)
            loc3 = int(primer_loc5) + len(primer) - 1 if strand == "0" \
                   else int(primer_loc5) - len(primer) + 1
            primer_info[primer] = "|".join(
                [readset, primer,strand,chrom,primer_loc5,str(loc3)])
            molecule = "-".join([chrom, strand, umi_loc, umi])

            if duplex_by_umi[molecule]["CC"] > 0:
                primer_metrics[primer]["CC_UMI"] += 1
                
            if duplex_by_umi[molecule]["TT"] > 0:
                primer_metrics[primer]["TT_UMI"] += 1
                
            if duplex_by_umi[molecule]["NN"] > 0:
                primer_metrics[primer]["NN_UMI"] += 1                
            
            if duplex_by_umi[molecule]["NN"] > 0 and \
               duplex_by_umi[molecule]["CC"] == 0 and \
               duplex_by_umi[molecule]["TT"] == 0: # UMIs with only NN tags
                primer_metrics[primer]["only_NN"] += 1
                
            if duplex_by_umi[molecule]["CC"] > 1 and duplex_by_umi[molecule]["TT"] > 1:
                primer_metrics[primer]["duplexUMI"] += 1
                
            primer_metrics[primer]["UMI"] += 1

            read_frags, remainder = divmod(duplex_by_umi[molecule]["CC"],2)
            assert remainder == 0, "Read accounting error !"
            primer_metrics[primer]["CC"] += read_frags
            
            read_frags, remainder = divmod(duplex_by_umi[molecule]["TT"],2)
            assert remainder == 0, "Read accounting error !"
            primer_metrics[primer]["TT"] += read_frags
            
            read_frags, remainder = divmod(duplex_by_umi[molecule]["NN"],2)
            assert remainder == 0, "Read accounting error !"            
            primer_metrics[primer]["NN"] += read_frags
            
            

    # write to output file
    with open(outfile, 'w') as OUT:
        OUT.write(out_header)
        OUT.write("\n")
        for primer in primers:    
            CC_umi     = primer_metrics[primer]["CC_UMI"]
            TT_umi     = primer_metrics[primer]["TT_UMI"]
            NN_umi     = primer_metrics[primer]["NN_UMI"]
            only_NN_umi =  primer_metrics[primer]["only_NN"]

            
            CC_read_frags  = primer_metrics[primer]["CC"]
            TT_read_frags = primer_metrics[primer]["TT"]
            NN_read_frags = primer_metrics[primer]["NN"]

            total_read_frags  = primer_metrics[primer]["read_frags"]            
            total_umi = primer_metrics[primer]["UMI"]
            duplex_umi = (CC_umi + TT_umi) - (total_umi - only_NN_umi)

            assert duplex_umi == primer_metrics[primer]["duplexUMI"]
            total_read_frags = CC_read_frags + TT_read_frags + NN_read_frags
            #assert total_read_frags == CC_read_frags + TT_read_frags + NN_read_frags, "Read accounting error !"
            
            out = [primer_info[primer], str(total_umi), str(total_read_frags), str(CC_umi), \
                   str(CC_read_frags), str(TT_umi), str(TT_read_frags), str(NN_umi), str(NN_read_frags), \
                   str(only_NN_umi), str(duplex_umi)]
            OUT.write("|".join(out))
            OUT.write("\n")


if __name__ == '__main__':
    readset   = sys.argv[1]
    paramfile = sys.argv[2]
    import core.run_config
    cfg = core.run_config.run(readset, paramfile)
    run(cfg)
