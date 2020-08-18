#
#	NGSPIPELINE - pipeline module
#
#	* Pipeline adapted from a paper in Nature
#	* comparing whole genome sequence variant
#	* callers. Full citation below. 
#	* 
#	* Hwang, K., Lee, I., Li, H. et al.
#	* Comparative analysis of whole-genome sequencing 
#	* pipelines to minimize false negative findings. 
#	* Sci Rep 9, 3219 (2019). 
#	* https://doi.org/10.1038/s41598-019-39108-2
#
#	Hwang, K., Lee, I., Li, H. et al.
#
###########################################################

#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh
source ${pip_dir}/modules/checkparams.mod.sh
source ${pip_dir}/modules/ngspipeline.mod.sh

workflow() {
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]} ["tmp_prefix"]=${argv[2]})

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	custom_call align_reads_to_reference "aligning reads to the reference with ${inputs['aligner_id']}..."

	custom_call call_variants "calling variants with ${inputs['caller_id']}..."
}

#
#  tasks
#
initialize_inputs_hash() {
	local status=0

	# 1. set default parameter values
	printf '  setting default parameter values...'
	# aligner() params:
	inputs['aligner_id']=NULL;			# name of he chosen aligner tool
	inputs['caller_id']=NULL;				# name of the chosen variant caller tool
	inputs['ref_base']=NULL;			# basename of the reference sequence file
	inputs['sample_id']=NULL;			# prefix for outputs
	inputs['fastq1']=NULL;				# input fastq file (1)
	inputs['fastq2']=NULL;				# input fastq file (2)
	inputs['threads']=1;				# number of threads to parallelize over
	inputs['maxmem']=100M;				# max memory per thread (approximately)
	inputs['recal_realign_on']='no';	# do recalibrations (options "yes", "no", "both")
	echo '...done.'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	value_from_json ${inputs['input_json']} '.aligner_id'	inputs['aligner_id']
	value_from_json ${inputs['input_json']} '.caller_id'	inputs['caller_id']
	value_from_json ${inputs["input_json"]} '.ref_base'		inputs["ref_base"]
	value_from_json ${inputs['input_json']} '.sample_id'	inputs['sample_id']
	value_from_json ${inputs['input_json']} '.fastq1_base'	inputs['fastq1_base']
	value_from_json ${inputs['input_json']} '.fastq2_base'	inputs['fastq2_base']
	value_from_json ${inputs['input_json']} '.threads'		inputs['threads']
	value_from_json ${inputs['input_json']} '.maxmem'		inputs['maxmem']
	value_from_json ${inputs['input_json']} '.recal_realign_on'		inputs['recal_realign_on']

	inputs['ref']=${ref_dir}/${inputs['ref_base']}
	inputs["fastq1"]=${rds_dir}/${inputs["fastq1_base"]}
	inputs["fastq2"]=${rds_dir}/${inputs["fastq2_base"]}
	inputs['prefix']="${inputs['sample_id']}.${inputs['aligner_id']}.${inputs['caller_id']}.${inputs['ref_base']%.fa}";	# prefix for output files
	inputs['log_file']="${log_dir}/${inputs['prefix']}.log";							# output lof file path 
	inputs['idx']="${ref_dir}/${inputs['aligner_id']}.${inputs['ref_base']}";				# aligner index file
	[[ $status == 0 ]] && echo '...done'

	# 3. check that inputs make sense
	printf '  checking that parameter values make sense...'

	## check that all required tools are installed
	### aligner

	[[ $status == 0 ]] && echo '...done'

	# 4. set up output logging and temporary files
	#set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }
}

#
#  run workflow
#
workflow "$@"