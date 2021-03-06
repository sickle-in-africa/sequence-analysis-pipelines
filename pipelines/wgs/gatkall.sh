#
#  gatkall
#	* NGS variant calling pipeline
#	* using the gatk 
#
#	Kevin Esoh
#
####################################
#!/usr/bin/env bash

source includes/locations.sh
source ${pro_dir}/includes/config.sh
source ${pro_dir}/includes/utilities.sh
source ${pip_dir}/wgs/modules/checkparams.mod.sh
source ${pip_dir}/wgs/modules/geneMapNGS.mod.sh


workflow() {
	local \
		argv=("$@") \
	
	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]} ['tmp_prefix']=${argv[2]})
	
	inputs['aligner_id']=bwa
	inputs['caller_id']=gatk

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	custom_call fq "checking read file quality..."

	custom_call trim "trimming read files..."

	custom_call gmap "performing Mapping/Alignment with GATKv4 and BWA ..."

	custom_call varcall "calling variants with gatk..."
}


#
#  run workflow
#
workflow "$@"
