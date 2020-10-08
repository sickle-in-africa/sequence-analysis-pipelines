#
#  SIMULATE COHORT - simulation pipeline
#
#	* generating a cohort of artificial reads 
#	* from an input reference sequence
#
#	Jack Morrice
#
###################################################
#!/bin/bash

source includes/locations.sh
source ${pro_dir}/includes/config.sh
source ${pro_dir}/includes/utilities.sh
source ${pip_dir}/wgs/modules/checkparams.mod.sh


workflow() { 
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} ["log_prefix"]=${argv[1]} ['tmp_prefix']=${argv[2]} )

	custom_call check_input_json "checking simulation input json file was provided..."

	custom_call initialize_inputs_hash "initializing simulation input parameter values..."

	custom_call simulate_cohort_reads "simulating reads for the cohort..."

	rm ${inputs['tmp_prefix']}*
}


#
#  tasks
#
initialize_inputs_hash() {
	local status=0

	# 1. set default parameter values
	printf '  setting default parameter values...'
	inputs["cohort_id"]="c_0"
	inputs["ref"]=NULL
	inputs["read_type"]="pe"; # pe: paried end, sr: single read
	inputs["n_samples"]=1
	inputs["threads"]=1
	echo '...done'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	_update_input_value_from_json 'cohort_id'
	_update_input_value_from_json 'ref_label'
	_update_input_value_from_json 'read_type'
	_update_input_value_from_json 'n_samples'
	_update_input_value_from_json 'threads'

	# 3. set derivative parameters
	inputs['ref']=${ref_dir}/${inputs['ref_label']}/${inputs['ref_label']}.fa.gz
	echo '...done'

	# 4. check that inputs make sense
	printf '  checking that parameter values make sense...'
	check_id "cohort" ${inputs["cohort_id"]} || { echo 'invalid choice of cohort id'; status=1; }
	check_ref || status=1
	check_int ${inputs["n_samples"]} n_samples || status=1
	[[ $status == 0 ]] && echo '...done'

	# 5. set up logging information
	set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }

	return $status
}

_update_input_value_from_json() {
	local key=$1
	value_from_json ${inputs['input_json']} '.'${key}	inputs[${key}]
}

simulate_cohort_reads() {

	custom_call_2 _mutate_reference '  mutating the reference...' || return 1

	custom_call_2 _generate_reads '  simulating reads...' || return 1

	custom_call_2 _sort_gvcf_file '  sorting truth gvcf files...' truth || return 1
}

_mutate_reference() {

	$bgzip -cd ${inputs["ref"]} > ${inputs['tmp_prefix']}${inputs['ref_label']}.fa

	$simulate MutateReference \
		--input_ref ${inputs['tmp_prefix']}${inputs['ref_label']}.fa \
		--output_ref ${inputs['tmp_prefix']}${inputs["cohort_id"]}.fa \
		--input_json ${inputs["input_json"]} \
		--output_vcf ${vcf_dir}/${inputs["cohort_id"]}.truth.g.vcf \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log

}

_generate_reads() {
	local \
		s \
		prefix \
		sample_log_file \
		option_string \
		log_file_string

	# 1. write sample id list file

	echo "# sample file for cohort: ${inputs["cohort_id"]}" > ${rds_dir}/${inputs["cohort_id"]}.samples.list
	for s in $(seq 1 1 ${inputs["n_samples"]}); do
		prefix="${inputs["cohort_id"]}.s_$s"

		# create a new log file for each sample if they don't already exist:
		sample_log_file=${inputs["log_prefix"]}${prefix}.log
		[[ -s $sample_log_file ]] || > $sample_log_file

		echo -e "s_$s\t${prefix}.raw_1P.fq.gz\t${prefix}.raw_2P.fq.gz\tILLUMINA" >> ${rds_dir}/${inputs["cohort_id"]}.samples.list
	done

	# 2. generate reads files

	option_string="GenerateReads \
		--input_ref ${inputs['tmp_prefix']}${inputs["cohort_id"]}.fa \
		--reads_prefix "${rds_dir}/"${inputs["cohort_id"]}.{1}.raw \
		--input_json ${inputs["input_json"]}"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{1}.log"

	run_in_parallel \
		$simulate \
		${rds_dir}/${inputs["cohort_id"]}.samples.list \
		"${option_string}" \
		${inputs["threads"]} \
		"${log_file_string}"

	# 3. zip reads files

	run_in_parallel \
		gzip \
		${rds_dir}/${inputs["cohort_id"]}.samples.list \
		"-f \
			${rds_dir}/${inputs['cohort_id']}.{1}.raw_1P.fq \
			${rds_dir}/${inputs['cohort_id']}.{1}.raw_2P.fq" \
		${inputs["threads"]} \
		"${log_file_string}" 
}

_sort_gvcf_file() {
	local input_gvcf_stage=$1

	$bcftools sort \
		${vcf_dir}/${inputs["cohort_id"]}.${input_gvcf_stage}.g.vcf \
		-o ${vcf_dir}/${inputs["cohort_id"]}.${input_gvcf_stage}.sorted.g.vcf \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log
}


#
#  run workflow
#
workflow "$@"
