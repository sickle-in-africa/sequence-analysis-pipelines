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

source includes/locations.sh
source ${pro_dir}/includes/utilities.sh
source ${pip_dir}/wgs/modules/checkparams.mod.sh
source ${pip_dir}/wgs/modules/ngspipeline-align.mod.sh
source ${pip_dir}/wgs/modules/ngspipeline-vcall.mod.sh

workflow() {
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]} ["tmp_prefix"]=${argv[2]})

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	custom_call align_reads_to_reference "aligning reads to the reference..."

	custom_call call_variants "calling variants..."

	custom_call clean_temp 'cleaning temporary files...'
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
	inputs['cohort_id']=NULL;			# prefix for outputs
	inputs['threads']=1;				# number of threads to parallelize over
	inputs['maxmem']=100M;				# max memory per thread (approximately)
	inputs['recal_realign_on']='no';	# do recalibrations (options "yes", "no", "both")
	echo '...done.'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	_update_input_value_from_json 'aligner_id'
	_update_input_value_from_json 'caller_id'
	_update_input_value_from_json 'ref_base'
	_update_input_value_from_json 'cohort_id'
	_update_input_value_from_json 'threads'
	_update_input_value_from_json 'maxmem'
	_update_input_value_from_json 'recal_realign_on'

	# 3. set derivative parameters
	inputs['ref']=${ref_dir}/${inputs['ref_base']}
	inputs['prefix']="${inputs['cohort_id']}.${inputs['aligner_id']}.${inputs['caller_id']}.${inputs['ref_base']%.fa}";	# prefix for output files
	inputs['log_file']="${log_dir}/${inputs['prefix']}.log";							# output lof file path 
	inputs['idx']="${ref_dir}/${inputs['aligner_id']}.${inputs['ref_base']}";			# aligner index file
	inputs['samples_list']=${rds_dir}/${inputs['cohort_id']}.samples.list
	[[ $status == 0 ]] && echo '...done'

	# 4. check that inputs make sense
	printf '  checking that parameter values make sense...'

	## check that all required tools are installed
	### aligner

	[[ $status == 0 ]] && echo '...done'

	# 5. set up output logging and temporary files
	set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }
}

_update_input_value_from_json() {
	local key=$1
	value_from_json ${inputs['input_json']} '.'${key}	inputs[${key}]
}

#
#  run workflow
#
workflow "$@"