#
#  BENCHMARK - test suite
#
#	* comparing multiple pipelines
#	* via comparison plots of speed
#	* and accuracy
#
#	Jack Morrice
#
#########################################
#!/bin/bash
cd '../'

source includes/locations.sh
source includes/utilities.sh

source ${tst_dir}/modules/test-suites.mod.sh
source ${pip_dir}/wgs/modules/checkparams.mod.sh

workflow() {
	local 
		argv=("$@")

	declare -A inputs=( ["input_json"]=${tst_dir}/${argv[0]} ["log_prefix"]=${argv[1]} ['tmp_prefix']=${argv[2]} )

	custom_call_ts check_input_json "checking a test suite input json file was provided..."

	custom_call_ts initialize_inputs_hash "initializing input parameter values..."

	custom_call_ts run_benchmark_test "running pipeline benchmark test..."

	custom_call_ts plot_benchmark_test "plotting benchmark test results..."

}


#
#  tasks
#
run_benchmark_test() {
	local \
		jacc_value='NULL'

	> ${dat_dir}/benchmark-runs.dat
	echo -e "run\tpipeline\tjaccard\truntime" >> ${dat_dir}/benchmark-runs.dat

	for i in $(seq 1 ${inputs["runs_per_pipeline"]}); do

		# clean data directories
		(cd ${tst_dir} && ./clean-data.sh ${inputs['cohort_id']}) || return 1

		custom_call_ts_2 simulate_reads "[run $i] simulating cohort reads..." || return 1

		for j in "gatkall" "bcfall"; do

			inputs["pipeline_id"]=$j
			inputs['aligner_id']='bwa'
			if [[ ${inputs["pipeline_id"]} == "gatkall" ]]; then
				inputs['caller_id']='gatk'
			elif [[ ${inputs['pipeline_id']} == 'bcfall' ]]; then
				inputs['caller_id']='samtools'
			fi
			inputs['prefix']="${inputs['cohort_id']}.${inputs['aligner_id']}.${inputs['caller_id']}.${inputs['ref_base']%.fa}";	# prefix for output files

			# clean data directories
			custom_call_ts_2 clean_downstream_data 'cleaning downstream pipeline data...' || return 1

			# 2. run the pipeline and time it:
			custom_call_ts_2 run_pipeline "[run $i] running pipeline: ${inputs['pipeline_id']} with aligner: ${inputs['aligner_id']} and caller: ${inputs['caller_id']}..." || return 1

			# 3. compare truth and estimated gvcf files

			custom_call_ts_2 compare_truth_est_vcf "[run $i] comparing the simulation's truth vcf file against the variant caller output..." || return 1
			
			custom_call_ts_2 jaccard_index "[run $i] computing jaccard index..." jacc_value || return 1
			
			echo "Jaccard value: $jacc_value"
			echo "pipeline runtime: $DIFF seconds"

			echo -e "$i\t${inputs["pipeline_id"]}\t${jacc_value}\t$DIFF" >> ${dat_dir}/benchmark-runs.dat

		done

	done
}

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
		|| { echo 'test suite ERROR: pipeline step failed'; exit 1; }
	END=$(date +%s.%N)
	DIFF=$(echo "$END - $START" | bc)
}

plot_benchmark_test() {

	Rscript $benchmarkR \
		${dat_dir}/benchmark-runs.dat \
		${tls_dir}/jaccard \
		${med_dir}/benchmark

}

#
#  run workfow
#
workflow "$@"
