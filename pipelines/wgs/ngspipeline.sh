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
source ${pro_dir}/includes/config.sh
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
	inputs['caller_id']=NULL;			# name of the chosen variant caller tool
	inputs['ref_label']=NULL;			# basename of the reference sequence file
	inputs['cohort_id']=NULL;			# prefix for outputs
	inputs['threads']=1;				# number of threads to parallelize over
	inputs['maxmem']=100M;				# max memory per thread (approximately)
	inputs['recal_realign_on']='no';	# do recalibrations (options "yes", "no", "both")
	echo '...done.'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	_update_input_value_from_json 'aligner_id'
	_update_input_value_from_json 'caller_id'
	_update_input_value_from_json 'ref_label'
	_update_input_value_from_json 'cohort_id'
	_update_input_value_from_json 'threads'
	_update_input_value_from_json 'maxmem'
	_update_input_value_from_json 'recal_realign_on'

	# 3. set derivative parameters
	inputs['ref']=${ref_dir}/${inputs['ref_label']}/${inputs['ref_label']}.fa.gz
	inputs['prefix']="${inputs['cohort_id']}.${inputs['aligner_id']}.${inputs['caller_id']}.${inputs['ref_label']}";	# prefix for output files
	inputs['log_file']="${log_dir}/${inputs['prefix']}.log";							# output lof file path 
	inputs['idx']="${ref_dir}/${inputs['aligner_id']}.${inputs['ref_label']}";			# aligner index file
	inputs['samples_list']=${rds_dir}/${inputs['cohort_id']}.samples.list
	[[ $status == 0 ]] && echo '...done'

	# 4. check that inputs make sense
	printf '  checking that parameter values make sense...'
	check_sample || status=1
	check_int ${inputs["threads"]} threads || status=1
	check_int ${inputs["trim_minlen"]} trim_minlen || status=1
	check_int ${inputs["trim_leadx"]} trim_leadx || status=1
	check_int ${inputs["trim_trailx"]} trim_trailx || status=1
	check_ref || status=1
	check_bwa_idx || status=1
	check_gatk_dict || status=1
	check_samtools_fai || status=1
	[[ $status == 0 ]] && echo '...done'

	printf '  checking that the necessary tools have been installed...'
	_check_aligner_installed || status=1
	_check_caller_installed  || status=1
	_check_samtools          || status=1
	_check_gatk              || status=1
	[[ $status == 0 ]] && echo '...done'

	# 5. set up output logging and temporary files
	set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }
}

_update_input_value_from_json() {
	local key=$1
	value_from_json ${inputs['input_json']} '.'${key}	inputs[${key}]
}

function check_sample() {
	local \
		line \
		line_array \
		sample_id \
		sample_log_file

	# check meta variabe was set
	if [[ "${inputs['samples_list']}" == NULL ]]; then
		echo -e "\e[38;5;1mERROR\e[0m: -s,--sample_list not provided! Exiting..."; 1>&2;
		return 1
	elif [ -f ${inputs['samples_list']} -a -s ${inputs['samples_list']} ]; then
		for i in $(awk '{print $1}' ${inputs['samples_list']}); do
			[ ! $i == "#"* ] && continue
			if [ ! -f ${rds_dir}/$i ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' was not found in the directory '${rds_dir}'.\nPlease specify the path with -p or --path or check that the files in the path are the same in the sample list" 1>&2;
				return 1;
			elif [ -f ${rds_dir}/${i} -a ! -s ${rds_dir}/${i} ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' may be empty. Please check and correct '${rds_dir}'." 1>&2;
				return 1;
		   fi
		done
	elif [ -f ${inputs['samples_list']} -a ! -s ${inputs['samples_list']} ]; then
		echo -e "\e[38;5;1mERROR\e[0m: '$meta' seems to be empty! Please check and correct." 1>&2;
	fi

	# check read files in $meta file exist, and create log files if they do not exist
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			sample_id=${line_array[0]}
			[[ -s ${rds_dir}/"${line_array[1]}" ]] || { echo "Error: ${line_array[0]} does not exist or is empty"; return 1; }
			[[ -s ${rds_dir}/"${line_array[2]}" ]] || { echo "Error: ${line_array[0]} does not exist or is empty"; return 1; }


		fi
	done < ${inputs['samples_list']}
}

_check_aligner_installed() {
	if [ ${inputs['aligner_id']} = "bwa" ]; then
		_check_bwa || { echo "error: bwa not installed"; return 1; }
	
	elif [ ${inputs['aligner_id']} = "gsnap" ]; then
		_check_gsnap || { echo "error: gsnap not installed"; return 1; }
	
	elif [ ${inputs['aligner_id']} = "bowtie2" ]; then
		_check_bowtie2 || { echo "error: bowtie2 not installed"; return 1; }
	
	elif [ "${inputs['aligner_id']}" = "stampy" ]; then
		_check_stampy || { echo "error: stampy not installed"; return 1; }
	fi
}

_check_caller_installed() {
	if [ "${inputs['caller_id']}" = "freebayes" ];then
		_check_freebayes || { echo "error: freebayes not installed"; return 1; }
	elif [ "${inputs['caller_id']}" = "samtools" ];then
		_check_bcftools || { echo "error: bcftools not installed"; return 1; }
	elif [ "${inputs['caller_id']}" = "gatk" ];then
		_check_gatk || { echo "error: gatk not installed"; return 1; }
	fi
}

_check_bowtie2() { $bowtie2 --help &> /dev/null || return 1; }
_check_samtools() {	$samtools --help &> /dev/null || return 1; }
_check_bcftools() {	$bcftools --help &> /dev/null || return 1; }
_check_bwa() { man $bwa.1 &> /dev/null || return 1; }
_check_gsnap() { $gsnap --help &> /dev/null || return 1; }
_check_stampy() { python $stampy --help &> /dev/null || return 1; }
_check_gatk() { $gatk --help &> /dev/null || return 1; }
_check_fastqc() { $fastqc --help &> /dev/null || return 1; }
_check_trimmomatic() { java -jar $trimmomatic -version &> /dev/null || return 1; }
_check_freebayes() { $freebayes --help &> /dev/null || return 1; }

#
#  run workflow
#
workflow "$@"
