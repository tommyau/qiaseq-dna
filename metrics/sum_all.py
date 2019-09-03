import os

def run(cfg):
    # get read set name
    readSet = cfg.readSet
    
    umi_filter = "umi_filter"
    duplex     = "duplex"
    prep       = "prep"
    if cfg.outputDetail: # more details if needed
        prep       = "prep.detail" if os.path.exists(readSet + ".prep.detail.summary.txt") else "prep" # ion-reads don't have prep.detail file
        umi_filter = "umi_filter.detail"
        duplex     = "duplex.detail"

    # concatenate summary files
    with open(readSet + ".sum_all.summary.txt", "w") as OUT:
        for fileType in (prep, "align", umi_filter, "umi_frags", "sum.uniformity.primer", "umi_depths", duplex, "smCounter", "vcf_complex"):
            '''
            prep, align           -->  trimming metrics                 : core/prep.py, misc/process_ion.py (for ion reads)
            umi_filter            -->  primer finding metrics           : core/umi_filter.py
            umi_frags             -->  reads per umi info               : metrics/umi_frags.py
            sum.uniformity.primer -->  primer level uniformity metrics  : metrics/sum_uniformity_primer.py
            umi_depths            -->  umi depth and lod info           : metrics/umi_depths.py
            duplex                -->  duplex level metrics             : metrics/duplex_summary.py
            smCounter             -->  variant calling threshold        : core/sm_counter_wrapper.py
            vcf_complex           -->  metrics from reconstructing 
                                       complex variants from primitives : annotate/vcf_complex.py
            '''
            fileName = readSet + "." + fileType + ".summary.txt"
            if os.path.isfile(fileName):
                for line in open(fileName):
                    OUT.write(line)
