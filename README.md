# qiaseq-dna
This repository contains some example code for processing reads from QIAGEN QIAseq DNA enrichment kits.

python run examples
-------------------
The "run" directory contains three examples of how to use the read processing code in this repository:

1. **run_dedup** - Use the Picard MarkDuplicates utility from the Broad Institute to remove PCR duplicate reads using the UMI and both paired-end read start locations on the reference genome.  This might be useful for subsequent germline variant calling, or structural variant detection applications.  Note that this script does not filter ligation chimera reads, nor does it remove PCR duplicate reads caused by internal re-priming by a downstream SPE primer. 

2. **run_consensus** - Use the "fgbio" package from Fulcrum Genomics to create a single consensus read pair for each input molecule.  This might be useful to prepare consensus read alignments for a subsequent SNP/indel variant calling procedure that does not use UMI-tagged reads (such as VarDict, MuTect, etc.).  Users might need to tune the fgbio parameters (in addtion to variant calling parameters) for their application.  We have not performed any variant calling performance benchmarking using the consensus-read BAM generated by this pipeline.

3. **run_sm_counter** - Use the "smCounter" variant calling procedure, described here:
[Detecting very low allele fraction variants using targeted DNA sequencing and a novel molecular barcode-aware variant caller", BMC Genomics, 2017 18:5](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3425-4)   
The smCounter variant caller uses a statistical model that requires both raw sequencer base calls and identification of unique input molecules using the UMI tag and the genome position of the random fragmentation site.

python packages
---------------------
#### core  
This package contains modules that trim common regions of the reads, align reads to the reference genome, identify putative original input molecules, and trim SPE primer regions from the genome alignments.  These steps generate a BAM file suitable for subsequent variant calling.

#### metrics
This package contains auxiliary modules that provide read accounting and enrichment summary metrics such as uniformity and fragment length distribution.

#### varcall
This package contains the smCounter variant caller and downstream VCF annotation using snpEff.

Docker image for third-party dependencies
-----------------------------------------
The python modules in this repository have many dependencies on third-party NGS software (e.g. BWA, samtools, etc.) and GNU Linux utilities (sort, zcat, etc).  Please **DO NOT ATTEMPT** to use the python modules in this git repository without first running the code on the example read set using our Docker image:

```bash
### Pull the docker image
sudo docker pull rpadmanabhan9/qiaseq-dna

### Run a container from the image above interactively, mounting your run directory i.e. directory where the output files will be created
sudo docker run -it -v /home/your_fav_dir/:/mnt/qiaseq-run/ rpadmanabhan9/qiaseq-dna

### Change directory and get the latest code from github
cd /srv/qgen/code/
git clone https://github.com/qiaseq/qiaseq-dna.git

### Change to run directory and copy over parameters file
cd /mnt/qiaseq-run/
cp /srv/qgen/code/qiaseq-dna/run_consensus.params.txt ./

### Edit the bottom of run_consensus.params.txt if you need to change the read set and primer file

### Run the pipeline
python /srv/qgen/code/qiaseq-dna/run_consensus.py  NEB_S2 run_consensus.params.txt  > run.log 2>&1 &  

```
The dependencies are fully documented in the Dockerfile in this repository.

Please address questions to raghavendra.padmanabhan@qiagen.com, with CC to john.dicarlo@qiagen.com.

