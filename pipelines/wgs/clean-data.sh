#  clean-data
#
#	* Clean all intermediate and final WGS PIPELINE output data
#	* for a cohort, *downstream of the read files*
#	* This script does not remove the raw read files. 
#
#	Jack Morrice
#
###############################################################

#!/bin/bash

source includes/locations.sh
source ${pro_dir}/includes/utilities.sh

workflow() { 
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]} ["tmp_prefix"]=${argv[2]})

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	if [[ -z ${inputs['cohort_id']} ]]; then
		echo 'Error: cohort id not supplied. Exiting...'
		exit 1
	fi

	echo "  cleaning ${inputs['cohort_id']} fastqc files" && rm ${fqc_dir}/${inputs['cohort_id']}* 2> /dev/null
	echo "  cleaning ${inputs['cohort_id']} trimmed read files" && rm ${rds_dir}/${inputs['cohort_id']}*.paired* ${rds_dir}/*.unpaired* 2> /dev/null
	echo "  cleaning ${inputs['cohort_id']} sam files" && rm ${sam_dir}/${inputs['cohort_id']}* 2> /dev/null
	echo "  cleaning ${inputs['cohort_id']} bam files" && rm ${bam_dir}/${inputs['cohort_id']}* 2> /dev/null
	echo "  cleaning ${inputs['cohort_id']} temporary files" && rm ${tmp_dir}/* 2> /dev/null

	echo "  cleaning ${inputs['cohort_id']} estimated vcf files"
	rm -r ${vcf_dir}/${inputs['cohort_id']}*.isec 2> /dev/null
	rm $(find ${vcf_dir} -maxdepth 1 -name ${inputs['cohort_id']}* | grep -v 'truth') 2> /dev/null

	#echo "cleaning log files" && rm -r ${log_dir}/* 2> /dev/null

	exit 0
}

#
#  tasks
#
initialize_inputs_hash() {

	# 1. set default parameter values
	printf '  setting default parameter values...'
	inputs["cohort_id"]=NULL
	echo '...done'

	printf '  updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.cohort_id'   inputs["cohort_id"]
	echo '...done.'

}

#
#  run workflow
#
workflow "$@"

