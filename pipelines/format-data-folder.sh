#
#  FORMAT DATA FOLDER - setup script for pipelines
#
#	* script to create all sub-directories within
#	* the data/ folder needed for the pipelines
#	* to run. 
#
#	Jack Morrice
#
##########################################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

workflow() {
	local argv=("$@")

	custom_call check_directories_exist "checking parent directories exist..."

	custom_call create_sub_directories	"creating data/ sub-directories..."
}


#
#  tasks
#
check_directories_exist() {
	[[ -d ${pro_dir} ]] || { echo 'specified project directory was not found. Please check includes/locations.sh and try again'; return 1; }
	[[ -d ${dat_dir} ]] || { echo 'specified data directory was not found. Please check includes/locations.sh and try again'; return 1; }
}

create_sub_directories() {
	local status=0

	[[ ! -d ${bam_dir} ]] &&  { mkdir ${bam_dir}; echo "  created ${bam_dir}"; }
	[[ ! -d ${fqc_dir} ]] &&  { mkdir ${fqc_dir}; echo "  created ${fqc_dir}"; }
	[[ ! -d ${log_dir} ]] &&  { mkdir ${log_dir}; echo "  created ${log_dir}"; }
	[[ ! -d ${med_dir} ]] &&  { mkdir ${med_dir}; echo "  created ${med_dir}"; }
	[[ ! -d ${rds_dir} ]] &&  { mkdir ${rds_dir}; echo "  created ${rds_dir}"; }
	[[ ! -d ${sam_dir} ]] &&  { mkdir ${sam_dir}; echo "  created ${sam_dir}"; }
	[[ ! -d ${tmp_dir} ]] &&  { mkdir ${tmp_dir}; echo "  created ${tmp_dir}"; }

	return $status
}

#
#  run workflow
#
workflow "$@"
