#
#  SAP - wrapper script
#		for sequence-analysis-pipelines library
#
#	* call all pipeline scripts with 
#	* this wrapper rather than the 
#	* pipeline scripts directly. 
#
#	* basic usage:
#	* ./sap.sh <study_type> <pipeline id> <pipeline inputs id>
#
#	* example:
#	* ./sap.sh wgs gatkall basic
#
##############################################################
#!/bin/bash

source includes/locations.sh
source includes/utilities.sh

workflow() {
	local argv=("$@")

	local \
		study_type=${argv[0]} \
		pipeline_id=${argv[1]}

	if [[ $# == 5 ]]; then
		local \
			pipeline_inputs_id=${argv[2]} \
			log_prefix=${argv[3]} \
			tmp_prefix=${argv[4]}

		${pip_dir}/${study_type}/${pipeline_id}.sh \
			${pin_dir}/${study_type}/${pipeline_id}.${pipeline_inputs_id}.json \
			$log_prefix \
			$tmp_prefix

	elif [[ $# == 3 ]]; then
		local pipeline_inputs_id=${argv[2]}

		${pip_dir}/${study_type}/${pipeline_id}.sh \
			${pin_dir}/${study_type}/${pipeline_id}.${pipeline_inputs_id}.json

	elif [[ $# == 2 ]]; then

		${pip_dir}/${study_type}/${pipeline_id}.sh

	else
		echo "sap error: invalid number of arguments supplied. basic usage is"
		echo "./sap.sh <study_type> <pipeline id> <pipeline inputs id>"
	fi

}

#
#  run workflow
#
workflow "$@"
