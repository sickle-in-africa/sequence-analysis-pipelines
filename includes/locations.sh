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
pro_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project
dat_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data
pip_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/pipelines
rds_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/reads
tmp_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/temp
fqc_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/fastqc
bam_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/bam
src_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/source
tst_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/test-suites
log_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/logs
tls_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools
med_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/media
ref_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/references
pin_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/inputs
vcf_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/vcf
sam_dir=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/data/sam

# tool executable paths
bwa=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/bwa/bwa
gatk=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/gatk-4.1.7.0/gatk
bcftools=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/bcftools-1.11/bin/bcftools
fastqc=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/FastQC/fastqc
generatejson=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/jaccard/generate-json.py
bgzip=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/htslib-1.11/bgzip
samtools=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/samtools-1.11/bin/samtools
benchmarkR=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/jaccard/benchmark.R
simulate=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/sequence-simulator/simulate.pl
jaccard=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/jaccard/jaccard.py
jq=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/jq-1.6/jq-linux64
trimmomatic=/mnt/lustre/groups/CBBI1243/SADaCC/sequence-analysis-projects/jmorrice/wgs-project/tools/Trimmomatic-0.39/trimmomatic-0.39.jar
