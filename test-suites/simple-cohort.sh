#
#	SIMPLE COHORT - bash pipeline
#
#	* compares truth and estimated vcf files
#	* for simulated SNP+indel data based on an input
#	* reference sequence, and outputs the Jaccard index
#	* for the truth and estimated vcf files.
#
#	Jack Morrice
#
################################################################
#!/bin/bash
cd '../'

source includes/locations.sh
source includes/utilities.sh

source ${tst_dir}/modules/test-suites.mod.sh
source ${pip_dir}/wgs/modules/checkparams.mod.sh

workflow () {
	local \
		argv=("$@") \
		jacc_value='NULL'

	declare -A inputs=( ["input_json"]=${tst_dir}/${argv[0]} ["log_prefix"]=${argv[1]} ['tmp_prefix']=${argv[2]} )

	custom_call_ts check_input_json "checking test suite input json file was provided..."

	custom_call_ts initialize_inputs_hash "initializing test suite input parameter values..."

	custom_call_ts simulate_reads 'simulating reads for a test cohort...'

	custom_call_ts run_pipeline 'running pipeline on simulated cohort read files...'

	custom_call_ts compare_truth_est_vcf "comparing the simulation's truth vcf file against the variant caller output..."
	
	custom_call_ts jaccard_index "computing jaccard index..." \
		jacc_value
	
	echo "  Jaccard value: $jacc_value"
	echo "  pipeline runtime: $DIFF seconds"

	custom_call_ts clean_temp "removing all temporary files..."
}

#
#  tasks
#
simulate_reads() {
	./sap.sh ${inputs['study_type']} ${inputs['simulate_id']} ${inputs['simulate_inputs_id']} \
		${inputs["log_prefix"]} \
		${inputs['tmp_prefix']} \
		|| return 1
}

run_pipeline() {
	START=$(date +%s.%N)
	./sap.sh ${inputs['study_type']} ${inputs["pipeline_id"]} ${inputs["pipeline_inputs_id"]} \
		${inputs["log_prefix"]} \
		${inputs['tmp_prefix']} \
		|| return 1
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)
}

#
#  run workflow
#
workflow "$@"