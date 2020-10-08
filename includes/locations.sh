#
#  LOCATIONS - bash config file
#
#	* file containing all path definitions
#	* for directories and executable tools
#	* used by the variant-caller package.
#	* 
#	* Use the <sap-path>/setup.sh script to 
#	* modify this file. Any direct changes
#	* to the file itself will be overwritten
#	* by another call of setup.sh.
#
#	Jack Morrice
#
################################################
	
# project path variables
pro_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test
dat_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data
pip_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/pipelines
rds_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/reads
tmp_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/temp
fqc_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/fastqc
bam_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/bam
src_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/source
tst_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/test-suites
log_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/logs
tls_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools
med_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/media
ref_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/references
pin_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/inputs
vcf_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/vcf
sam_dir=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/data/sam

# tool executable paths
bwa=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/bwa/bwa
gatk=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/gatk-4.1.7.0/gatk
bcftools=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/bcftools-1.11/bin/bcftools
fastqc=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/FastQC/fastqc
generatejson=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/jaccard/generate-json.py
bgzip=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/htslib-1.11/bgzip
samtools=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/samtools-1.11/bin/samtools
benchmarkR=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/jaccard/benchmark.R
simulate=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/sequence-simulator/simulate.pl
jaccard=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/jaccard/jaccard.py
jq=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/jq-1.6/jq-linux64
trimmomatic=/mnt/lustre/users/jmorrice/SADaCC/sequence-analysis-projects/development/wgs-test/tools/Trimmomatic-0.39/trimmomatic-0.39.jar
