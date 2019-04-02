from collections import defaultdict
import sys

import pysam
import numpy as np

# metrics
TOTAL_UMIS = 0
NUM_UMIS_NN_ONLY  = 1
NUM_READS_NN = 2
PERC_OF_UMIS_WITH_1_READ_FRAG = 3
PERC_OF_UMIS_WITH_2_READ_FRAGS = 4
PERC_OF_UMIS_WITH_3_READ_FRAGS = 5
NUM_UMIS_ALL_CC_1_READ_FRAG = 6
NUM_UMIS_ALL_TT_1_READ_FRAG = 7
NUM_UMIS_BOTH_CC_AND_TT_1_READ_FRAG_EACH = 8
NUM_UMIS_WITH_LT_2_CC_AND_TT_READ_FRAGS = 9

NUM_UMIS_GTE_2_READ_FRAGS_CC_LTE_1_READ_FRAG_TT = 10
NUM_UMIS_GTE_2_READ_FRAGS_TT_LTE_1_READ_FRAG_CC = 11

NUM_UMIS_GTE_2_READ_FRAGS_CC_0_READ_FRAG_TT = 12
NUM_UMIS_GTE_2_READ_FRAGS_TT_0_READ_FRAG_CC = 13

NUM_UMIS_ATLEAST_3_CC_OR_TT_READ_FRAGS = 14
NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS = 15

DUPLEX_UMIS = 16
DUPLEX_RATE = 17

MEAN_USABLE_UMI_RPU = 18
PERCENTILE_USABLE_UMI_RPU_25 = 19
PERCENTILE_USABLE_UMI_RPU_50 = 20
PERCENTILE_USABLE_UMI_RPU_75 = 21
MEAN_DUPLEX_UMI_RPU = 22
PERCENTILE_DUPLEX_UMI_RPU_25 = 23
PERCENTILE_DUPLEX_UMI_RPU_50 = 24
PERCENTILE_DUPLEX_UMI_RPU_75 = 25

NUM_METRICS_TOTAL = 26

def run(cfg):
    '''
    '''
    print('duplex_summary starting...')
    
    readset = cfg.readSet
    inbam   = cfg.readSet + '.umi_merge.bam'
   
    out_summary_file = readset + '.duplex.summary.txt'
    out_detail_file  = readset + '.duplex.detail.summary.txt'
    
    metric_vals     = [0]*NUM_METRICS_TOTAL

    IN = pysam.AlignmentFile(inbam,"rb")

    duplex_by_umi = defaultdict(lambda:defaultdict(int))
    umi_reads     = defaultdict(lambda:defaultdict(int))

    usable_umi_rpu = []
    duplex_umi_rpu = []

    for read in IN:
        umi =  read.get_tag(cfg.tagNameUmi)
        duplex_id = read.get_tag(cfg.tagNameDuplex)
        duplex_by_umi[umi][duplex_id] += 1


    # Calculate UMI counts by read fragments
    umis_with_1_readfrag = 0
    umis_with_2_readfrag = 0
    umis_with_3_readfrag = 0    
    for umi in duplex_by_umi:
        if duplex_by_umi[umi]['CC'] == 0 and duplex_by_umi[umi]['TT'] == 0 and duplex_by_umi[umi]['NN'] > 0:
            metric_vals[NUM_UMIS_NN_ONLY] += 1 # count NN only UMIs and skip
            continue
        
        total_reads_for_this_umi = 0
        for duplex, num_reads in duplex_by_umi[umi].iteritems():
            if duplex == 'NN':
                metric_vals[NUM_READS_NN] += 1
                continue # skip NN reads
            total_reads_for_this_umi += num_reads
            
        if total_reads_for_this_umi == 2: # 1 read frag
            umis_with_1_readfrag += 1
        elif total_reads_for_this_umi == 4: # 2 read frag
            umis_with_2_readfrag += 1
        elif total_reads_for_this_umi == 6: # 3 read frag
            umis_with_3_readfrag += 1

        total_reads_for_this_umi = 0
        
        assert total_reads_for_this_umi % 2 == 0, "Read accounting error !" # check to make sure we have 2 reads for each fragment
        

    for umi in duplex_by_umi:
        metric_vals[TOTAL_UMIS] += 1
        
        # check for NN only UMIs
        if duplex_by_umi[umi]['CC'] == 0 and duplex_by_umi[umi]['TT'] == 0 and duplex_by_umi[umi]['NN'] > 0:
            continue # skip
        
        if duplex_by_umi[umi]['CC'] <= 2 and duplex_by_umi[umi]['TT'] <= 2: # < 2 read frags for both CC and TT
            ratio_CC = float(duplex_by_umi[umi]['CC'])/(duplex_by_umi[umi]['TT']+duplex_by_umi[umi]['CC'])
            if ratio_CC == 1.0:
                metric_vals[NUM_UMIS_ALL_CC_1_READ_FRAG] +=1  # Singleton CC UMIs (1 read frag)
            elif ratio_CC == 0.0:
                metric_vals[NUM_UMIS_ALL_TT_1_READ_FRAG] += 1 # Singleton TT UMIs (1 read frag)
            else:
                metric_vals[NUM_UMIS_BOTH_CC_AND_TT_1_READ_FRAG_EACH] += 1 # UMIs with 1 CC and 1 TT read frag each
            metric_vals[NUM_UMIS_WITH_LT_2_CC_AND_TT_READ_FRAGS] += 1
            continue
        
        if duplex_by_umi[umi]['CC'] >= 4 or duplex_by_umi[umi]['TT'] >= 4:
            usable_umi_rpu.append(duplex_by_umi[umi]['CC']/2 + duplex_by_umi[umi]['TT']/2) # divide by 2 for read fragments
            metric_vals[NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS] += 1
        
        if duplex_by_umi[umi]['CC'] >= 6 or duplex_by_umi[umi]['TT'] >= 6:
            metric_vals[NUM_UMIS_ATLEAST_3_CC_OR_TT_READ_FRAGS] += 1
            
        if duplex_by_umi[umi]['CC'] >= 4 and duplex_by_umi[umi]['TT'] >= 4:
            duplex_umi_rpu.append(duplex_by_umi[umi]['CC']/2 + duplex_by_umi[umi]['TT']/2) # divide by 2 for read fragments
            metric_vals[DUPLEX_UMIS] += 1 # UMIs with CC >= 2 read frag and TT >= 2 read frag
            
        num_CC = 0 if duplex_by_umi[umi]['CC'] < 4 else duplex_by_umi[umi]['CC'] # treat as if no CC if < 2 read frag support for duplex
        num_TT = 0 if duplex_by_umi[umi]['TT'] < 4  else duplex_by_umi[umi]['TT'] # treat as if no TT if < 2 read frag support for duplex
        ratio_CC = float(num_CC)/float(num_CC + num_TT)
        if ratio_CC == 1.0:
            metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_CC_LTE_1_READ_FRAG_TT] += 1  # UMIs with TT = 1 read frag, CC >= 2 read frag
        elif ratio_CC == 0.0:
            metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_TT_LTE_1_READ_FRAG_CC]  += 1  # UMIs with CC = 1 read frag, TT >= 2 read frag

        num_CC = 0 if duplex_by_umi[umi]['CC'] == 0 else duplex_by_umi[umi]['CC'] # no CC, but >=2 read frag
        num_TT = 0 if duplex_by_umi[umi]['TT'] == 0 else duplex_by_umi[umi]['TT'] # no TT, but >=2 read frag
        ratio_CC = float(num_CC)/float(num_CC + num_TT)
        if ratio_CC == 1.0:
            metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_CC_0_READ_FRAG_TT] += 1  # UMIs with TT = 0 read frag, CC >= 2 read frag
        elif ratio_CC == 0.0:
            metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_TT_0_READ_FRAG_CC]  += 1  # UMIs with CC = 0 read frag, TT >= 2 read frag

   
    assert metric_vals[TOTAL_UMIS] == metric_vals[NUM_UMIS_ALL_CC_1_READ_FRAG] + \
        metric_vals[NUM_UMIS_ALL_TT_1_READ_FRAG] + metric_vals[NUM_UMIS_BOTH_CC_AND_TT_1_READ_FRAG_EACH] + \
        metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_CC_LTE_1_READ_FRAG_TT] + metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_TT_LTE_1_READ_FRAG_CC] + \
        metric_vals[DUPLEX_UMIS] + metric_vals[NUM_UMIS_NN_ONLY], \
        "Read accounting error !"

    assert metric_vals[TOTAL_UMIS] == metric_vals[NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS] + \
        metric_vals[NUM_UMIS_WITH_LT_2_CC_AND_TT_READ_FRAGS] + \
        metric_vals[NUM_UMIS_NN_ONLY], \
        "Read accounting error !"

    assert metric_vals[DUPLEX_UMIS] == metric_vals[NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS] - \
        metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_CC_LTE_1_READ_FRAG_TT] - metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_TT_LTE_1_READ_FRAG_CC], \
        "Read accounting error !"
    
    assert metric_vals[DUPLEX_UMIS] > 0, "Zero Duplex UMIs !"

    assert metric_vals[NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS] == len(usable_umi_rpu), "Read accounting error !"
    assert metric_vals[DUPLEX_UMIS] == len(duplex_umi_rpu), "Read accounting error !"

    metric_vals[DUPLEX_RATE] =  round(float(metric_vals[DUPLEX_UMIS]) / metric_vals[TOTAL_UMIS],2)
    # round some other values before printing
    metric_vals[PERC_OF_UMIS_WITH_1_READ_FRAG] = round(float(umis_with_1_readfrag)/metric_vals[TOTAL_UMIS],2)*100
    metric_vals[PERC_OF_UMIS_WITH_2_READ_FRAGS] = round(float(umis_with_2_readfrag)/metric_vals[TOTAL_UMIS],2)*100
    metric_vals[PERC_OF_UMIS_WITH_3_READ_FRAGS] = round(float(umis_with_3_readfrag)/metric_vals[TOTAL_UMIS],2)*100

    # get percentiles read frags per UMI for usable and duplex UMIs
    usable_umi_rpu  = np.array(usable_umi_rpu)
    duplex_umi_rpu  = np.array(duplex_umi_rpu)
    metric_vals[MEAN_USABLE_UMI_RPU] = round(np.mean(usable_umi_rpu),2)
    metric_vals[PERCENTILE_USABLE_UMI_RPU_25] = round(np.percentile(usable_umi_rpu,25,overwrite_input=True),2) # overwrite_input : array will get sorted only the first time
    metric_vals[PERCENTILE_USABLE_UMI_RPU_50] = round(np.percentile(usable_umi_rpu,50),2)
    metric_vals[PERCENTILE_USABLE_UMI_RPU_75] = round(np.percentile(usable_umi_rpu,75),2)
    metric_vals[MEAN_DUPLEX_UMI_RPU] = round(np.mean(duplex_umi_rpu),2)
    metric_vals[PERCENTILE_DUPLEX_UMI_RPU_25] = round(np.percentile(duplex_umi_rpu,25,overwrite_input=True),2)
    metric_vals[PERCENTILE_DUPLEX_UMI_RPU_50] = round(np.percentile(duplex_umi_rpu,50),2)
    metric_vals[PERCENTILE_DUPLEX_UMI_RPU_75] = round(np.percentile(duplex_umi_rpu,75),2)

    all_metrics = [
        (str(metric_vals[NUM_UMIS_NN_ONLY]), "UMIs with only NN"),
        (str(metric_vals[NUM_READS_NN]), "read fragments with NN excluded in analysis metrics"),
        (str(metric_vals[PERC_OF_UMIS_WITH_1_READ_FRAG]), "% of UMIs with 1 read frag "),
        (str(metric_vals[PERC_OF_UMIS_WITH_2_READ_FRAGS]), "% of UMIs with 2 read frag"),
        (str(metric_vals[PERC_OF_UMIS_WITH_3_READ_FRAGS]), "% of UMIs with 3 read frag"),
        (str(metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_CC_0_READ_FRAG_TT]), "UMIs with >= 2 read frags CC and 0 read frag TT"),
        (str(metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_TT_0_READ_FRAG_CC]), "UMIs with >= 2 read frags TT and 0 read frag CC"),
        (str(metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_CC_LTE_1_READ_FRAG_TT]), "UMIs with >= 2 read frags CC and <= 1 read frag TT"),
        (str(metric_vals[NUM_UMIS_GTE_2_READ_FRAGS_TT_LTE_1_READ_FRAG_CC]), "UMIs with >= 2 read frags TT and <= 1 read frag CC"),
        (str(metric_vals[NUM_UMIS_ATLEAST_3_CC_OR_TT_READ_FRAGS]), "UMIs with >= 3 read frags CC OR TT"),
        (str(metric_vals[NUM_UMIS_WITH_LT_2_CC_AND_TT_READ_FRAGS]), "UMIs with < 2 read frags CC AND TT"),
        (str(metric_vals[DUPLEX_UMIS]), "usable Duplex UMIs (>= 2 read frags CC and TT)"),
        (str(metric_vals[NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS]), "usable UMIs (>= 2 read frags CC or TT)"),        
        (str(metric_vals[DUPLEX_RATE]), "Duplex Rate (Duplex UMIs/Total UMI)"),
        (str(metric_vals[NUM_UMIS_ALL_TT_1_READ_FRAG]), "UMIs with all TT (1 read frag UMIs)"),
        (str(metric_vals[NUM_UMIS_ALL_CC_1_READ_FRAG]), "UMIs with all CC (1 read frag UMIs)"),
        (str(metric_vals[NUM_UMIS_BOTH_CC_AND_TT_1_READ_FRAG_EACH]),"UMIs with both CC and TT (1 read frag CC and TT)"),
        (str(metric_vals[MEAN_USABLE_UMI_RPU]),"read fragments per usable UMI, mean"),
        (str(metric_vals[PERCENTILE_USABLE_UMI_RPU_25]),"read fragments per usable UMI, 25th percentile"),
        (str(metric_vals[PERCENTILE_USABLE_UMI_RPU_50]),"read fragments per usable UMI, 50th percentile"),
        (str(metric_vals[PERCENTILE_USABLE_UMI_RPU_75]),"read fragments per usable UMI, 75th percentile"),
        (str(metric_vals[MEAN_DUPLEX_UMI_RPU]),"read fragments per usable Duplex UMI, mean"),
        (str(metric_vals[PERCENTILE_DUPLEX_UMI_RPU_25]),"read fragments per usable Duplex UMI, 25th percentile"),
        (str(metric_vals[PERCENTILE_DUPLEX_UMI_RPU_50]),"read fragments per usable Duplex UMI, 50th percentile"),
        (str(metric_vals[PERCENTILE_DUPLEX_UMI_RPU_75]),"read fragments per usable Duplex UMI, 75th percentile")
    ]    
    summary_metrics = [
        (str(metric_vals[DUPLEX_UMIS]), "usable Duplex UMIs (>= 2 read frags CC and TT)"),
        (str(metric_vals[NUM_UMIS_ATLEAST_2_CC_OR_TT_READ_FRAGS]), "usable UMIs (>= 2 read frags CC or TT)"),
        (str(metric_vals[DUPLEX_RATE]), "Duplex Rate (Duplex UMIs/UMIs)"),
    ]
    assert len(all_metrics) == NUM_METRICS_TOTAL - 1 # no total UMIs in metric output

    with open(out_summary_file, 'w') as OUT:
        for element in summary_metrics:
            OUT.write("\t".join(element))
            OUT.write("\n")
            
    with open(out_detail_file, 'w') as OUT:
        for e in all_metrics:
            row = "\t".join(e)
            OUT.write(row)
            OUT.write("\n")

    IN.close()


if __name__ == '__main__':
    readset   = sys.argv[1]
    paramfile = sys.argv[2]
    import core.run_config
    cfg = core.run_config.run(readset, paramfile)
    run(cfg)

