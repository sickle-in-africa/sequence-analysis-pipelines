#
#  LOCATIONS - bash config file
#
#	* file containing all path definitions
#	* for directories and executable tools
#	* used by the variant-caller package.
#	* 
#	* Use the <sap-path>/setup.sh script to modify this
#	* file. Any direct changes to the file 
#	* itself will be overwritten by another call of setup.sh.
#
#	Jack Morrice
#
################################################
	
# project path variables
ref_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/references
pro_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines
log_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/logs
sam_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/sam
pip_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/pipelines
med_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/media
fqc_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/fastqc
vcf_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/vcf
rds_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/reads
tls_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools
src_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/source
bam_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/bam
tmp_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data/temp
tst_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/test-suites
pin_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/inputs
dat_dir=/home/jack/Local/GeneMap/sequence-analysis-pipelines/data

# tool executable paths
trimmomatic=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/trimmomatic-0.39/trimmomatic-0.39.jar
gmap_build=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/gmap-2019-09-12/bin/gmap_build
gsnap=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/gmap-2019-09-12/bin/gsnap
benchmarkR=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/jaccard/benchmark.R
jaccard=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/jaccard/jaccard.py
samtools=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/samtools-1.10/samtools
freebayes=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/freebayes/bin/freebayes
stampy=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/stampy-1.0.32/stampy.py
simulate=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/simulate-0.1/simulate.pl
bcftools=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/bcftools-1.10.2/bcftools
gatk=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/gatk-4.1.7.0/gatk
bowtie2_build=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/bowtie2-2.4.1-linux-x86_64/bowtie2-build
fastqc=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/fastqc/fastqc
bowtie2=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/bowtie2-2.4.1-linux-x86_64/bowtie2
generatejson=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/jaccard/generate-json.py
bwa=/home/jack/Local/GeneMap/sequence-analysis-pipelines/tools/bwa-0.7.17/bwa
