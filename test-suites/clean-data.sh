#
#  CLEAN DATA - maintenence script for pipeline use
#
#	* clean all intermediate and final TEST SUITE output data
#
#	Jack Morrice
#
###############################################################

#!/bin/bash
cd '../'

source includes/locations.sh
source includes/utilities.sh

workflow() { 
	local argv=("$@")

	declare -A inputs=( ['cohort_id']=${argv[0]} ['study_type']=wgs)

	custom_call_ts_2 clean_downstream_data 'cleaning downstream pipeline data...' || return 1

	[[ -d ${rds_dir} ]] && echo "  cleaning ${inputs['cohort_id']} reads files" && rm ${rds_dir}/${inputs['cohort_id']}* 2> /dev/null
	[[ -d ${rds_dir} ]] && echo "  cleaning ${inputs['cohort_id']} truth vcf files" && rm ${vcf_dir}/${inputs['cohort_id']}*.truth.* 2> /dev/null
}

#
#  tasks
#
clean_downstream_data() {
	local \
		temp_json \
		temp_json_id

	temp_json_id=$(random_id)
	temp_json=${pin_dir}/${inputs['study_type']}/clean-data.${temp_json_id}.json

	echo '{ "cohort_id": "'${inputs['cohort_id']}'"}' > $temp_json

	./sap.sh wgs clean-data $temp_json_id || return 1
	rm $temp_json
}

#
#  run workflow
#
workflow "$@"