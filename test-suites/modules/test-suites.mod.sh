#
#	TEST SUITES - bash module
#
#	* module of bash functions for use in 
#	* constructing generic pipeline test
#	* suites. 
#
#	Jack Morrice
#
###########################################

#  summary:
#	1. initialize_inputs_hash()
#	2. ...

initialize_inputs_hash() {
	local \
		cohort_id_pi \
		cohort_id_si \
		n_samples \
		status=0

	# 1. set default parameter values
	printf '  setting default parameter values...'
	inputs["runs_per_pipeline"]=1
	inputs['study_type']=NULL
	inputs["cohort_id"]=NULL
	inputs["simulate_id"]=simulate-cohort
	inputs["simulate_inputs_id"]=basic
	inputs["pipeline_id"]=gatkall
	inputs["pipeline_inputs_id"]=basic
	echo '...done'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	_update_input_value_from_json 'runs_per_pipeline'
	_update_input_value_from_json 'study_type'
	_update_input_value_from_json 'cohort_id'
	_update_input_value_from_json 'simulate_id'
	_update_input_value_from_json 'simulate_inputs_id'
	_update_input_value_from_json 'pipeline_id'
	_update_input_value_from_json 'pipeline_inputs_id'
	echo '...done'

	# 3. set derivative parameters
	inputs['aligner_id']=NULL
	inputs['caller_id']=NULL
	inputs['ref_base']=NULL
	if [[ ${inputs['pipeline_id']} == 'gatkall' ]]; then
		inputs['aligner_id']='bwa'
		inputs['caller_id']='gatk'
	elif [[ ${inputs['pipeline_id']} == 'bcfall' ]]; then
		inputs['aligner_id']='bwa'
		inputs['caller_id']='samtools'
	elif [[ ${inputs['pipeline_id']} == 'ngspipeline' ]]; then
		value_from_json ${pin_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json '.aligner_id'	inputs['aligner_id']
		value_from_json ${pin_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json '.caller_id'	inputs['caller_id']
	fi

	value_from_json ${pin_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json '.ref_base'	inputs['ref_base']
	inputs['prefix']="${inputs['cohort_id']}.${inputs['aligner_id']}.${inputs['caller_id']}.${inputs['ref_base']%.fa}";	# prefix for output files

	# 4. check that inputs make sense
	printf '  checking that parameter values make sense...'
	check_int ${inputs["runs_per_pipeline"]} runs_per_pipeline || status=1
	check_id "cohort" ${inputs["cohort_id"]} || status=1
	## are all the input cohort ids the same?
	value_from_json \
		${pin_dir}/${inputs['study_type']}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json \
		'.cohort_id' \
		cohort_id_si
	value_from_json ${pin_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json \
		'.cohort_id' \
		cohort_id_pi
	[ "${inputs["cohort_id"]}" == "$cohort_id_si" -a "$cohort_id_si" == "$cohort_id_pi" ] \
		|| { echo 'test suite ERROR: cohort ids do not match, please check all input json files and try again'; status=1; }
	[ -s ${pip_dir}/${inputs['study_type']}/${inputs["simulate_id"]}.sh ] || { echo "test suite ERROR: ${inputs["simulate_id"]}.sh file not found or is empty"; status=1; }
	[ -s ${pin_dir}/${inputs['study_type']}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json ] || { echo "test suite ERROR: ${pip_dir}/${inputs['study_type']}/${inputs["simulate_id"]}.${inputs["simulate_inputs_id"]}.json file not found or is empty"; status=1; }
	[ -s ${pip_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.sh ] || { echo "test suite ERROR: ${inputs["pipeline_id"]}.sh file not found or is empty"; status=1; }
	[ -s ${pin_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json ] || { echo "test suite ERROR: ${pip_dir}/${inputs['study_type']}/${inputs["pipeline_id"]}.${inputs["pipeline_inputs_id"]}.json file not found or is empty"; status=1; }

	## still need to write checks for the other input parameters...
	[[ $status == 0 ]] && echo '...done'

	# 5. set up output logging and temporary files
	set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }

	return $status
}

_update_input_value_from_json() {
	local key=$1
	value_from_json ${inputs['input_json']} '.'${key}	inputs[${key}]
}

compare_truth_est_vcf()
{
	# compare truth and estimated vcf files
	# truth:     vcf file produced by the simulator
	# estimated: vcf file produced by the variant caller

	# need to zip the vcf files, because bcf tools only accepts
	# gzipped inputs

	printf "  zipping genotyped vcf files..."
	$bcftools view \
		${vcf_dir}/${inputs['prefix']}.genotyped.g.vcf \
		-Oz \
		-o ${vcf_dir}/${inputs['prefix']}.genotyped.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	printf "  zipping truth vcf files..."
	$bcftools view \
		${vcf_dir}/${inputs['cohort_id']}.truth.sorted.g.vcf \
		-Oz \
		-o ${vcf_dir}/${inputs['cohort_id']}.truth.sorted.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	printf "  building indices for the zipped vcf files..."
	$bcftools index ${vcf_dir}/${inputs['prefix']}.genotyped.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; }
	$bcftools index ${vcf_dir}/${inputs['cohort_id']}.truth.sorted.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	printf "  computing the set difference and intersection vcf files..."
	$bcftools isec \
		${vcf_dir}/${inputs['prefix']}.genotyped.g.vcf.gz \
		${vcf_dir}/${inputs['cohort_id']}.truth.sorted.g.vcf.gz \
		-p ${vcf_dir}/${inputs['prefix']}.isec \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }
	
}

jaccard_index() {
	local \
		tmp_prefix=${tmp_dir}/$(random_id)_ \
		__resultvar=$1 \
		myresult='NULL'

	# compute the jaccard index from the output of
	# bcftools isec
	# assumes:
	#	0000, 0001 are set difference files, and
	#	0002, 0003 are intersection files

	printf '  computing jaccard index...'
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs['prefix']}.isec/0000.vcf 1> ${inputs['tmp_prefix']}${inputs["cohort_id"]}.isec \
		2>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; }
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs['prefix']}.isec/0001.vcf 1>> ${inputs['tmp_prefix']}${inputs["cohort_id"]}.isec \
		2>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; }
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs['prefix']}.isec/0002.vcf 1>> ${inputs['tmp_prefix']}${inputs["cohort_id"]}.isec \
		2>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; }
	$gatk CountVariants \
		-V ${vcf_dir}/${inputs['prefix']}.isec/0003.vcf 1>> ${inputs['tmp_prefix']}${inputs["cohort_id"]}.isec \
		2>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; }

	myresult=$(python $jaccard \
		--input ${inputs['tmp_prefix']}${inputs["cohort_id"]}.isec)

	eval $__resultvar="'$myresult'"

	echo '...done'
}