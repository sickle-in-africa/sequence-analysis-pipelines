#
#  UTILITIES - bash module
#
#	* module of generic bash functions
#	* that have uses in scripts across the
#	* whole variant-caller package. 
#
#	Jack Morrice
#
##############################################

## output formatting
red='\033[0;31m'
green='\033[0;32m'
blue='\033[0;34m'
nc='\033[0m' # No Color
error="\e[38;5;1mERROR\e[0m:"

# generic useful functions

custom_call() {
	# for *pipeline top level functions*
	# call tasks above with colored output
	# and terminate the workflow on errors
	local \
		routine=$1 \
		msg=$2
	shift; shift

	printf "${green}$msg${nc}"; echo
 
	$routine "$@" \
		|| { printf "${red}...failed!${nc}"; echo; exit 1; } \
		&& { printf "${green}...done."${nc}; echo; }
}

custom_call_2() {
	# for *pipeline 2nd level functions*
	# call tasks above with messages
	# and return 1 on errors
	local \
		routine=$1 \
		msg=$2
	shift; shift
 
 	printf "$msg"
	if $routine "$@" ; then
		echo "...done"
	else
		echo "...failed!"
		return 1
	fi
}

custom_call_ts() {
	# for *test-suite top level functions*
	# call tasks above with colored output
	# and terminate the workflow on errors
	local \
		routine=$1 \
		msg=$2
	shift; shift

	printf "${blue}xxx $msg${nc}"; echo
 
	$routine "$@" \
		|| { printf "${red}xxx ...failed!${nc}"; echo; exit 1; } \
		&& { printf "${blue}xxx ...done."${nc}; echo; }
}

custom_call_ts_2() {
	# for *test-suite 2nd level functions*
	# call tasks above with colored output
	# and terminate the workflow on errors
	local \
		routine=$1 \
		msg=$2
	shift; shift

	printf "${blue}xxx+++ $msg${nc}"; echo
 
	$routine "$@" \
		|| { printf "${red}xxx+++ ...failed!${nc}"; echo; return 1; } \
		&& { printf "${blue}xxx+++ ...done."${nc}; echo; }
}

check_input_json() {

    # check a filename was given
    [[ -s ${inputs["input_json"]} ]] || { echo 'ERROR: please specify an input json file'; return 1; }

    # check that it is a valid json file
    cat ${inputs["input_json"]} | $jq -e . 1> /dev/null || { echo 'ERROR in input json file'; return 1; }
}

value_from_json() {
	local file=$1		# path of input json file
	local key=$2		# data field/key of {key, value} in the file
	local dest=$3		# destination to save value to
	
	# save the value, removing any quotation marks if present
	local value=$(sed -e 's/^"//' -e 's/"$//' <<<"$(cat $file | $jq $key)")

	# if the value is not null then the key exists in the json file, 
	# and we can update the $dest variable with $value
	[[ $value != "null" ]] && eval $dest=$value
}

random_id() {
	openssl rand -hex 4
}

check_int() {
	local value=$1
	local name=$2

	[[ $1 == 'NULL' ]] || return 0

	re=^[0-9]+$
	if [[ ! $1 =~ ^[0-9]+$ ]] ; then
   		echo "error: $2 is not a positive integer" >&2; return 1
	fi
}

run_in_parallel() {
	local \
		routine=$1 \
		parameter_file=$2 \
		option_string=$3 \
		n_threads=$4 \
		output_string=$5 

	sed '/^#/d' ${parameter_file} | \
		parallel \
			--col-sep '\t' \
			echo -e "${option_string} \&\>\> ${output_string}" | \
			xargs \
				-I input \
				-P $n_threads \
				bash -c "$routine input" \
				|| return 1
}

run_gatk_in_parallel() {
	local \
		java_options=$1 \
		parameter_file=$2 \
		option_string=$3 \
		n_threads=$4 \
		output_string=$5 

	sed '/^#/d' ${parameter_file} | \
		parallel \
			--col-sep '\t' \
			echo -e "${option_string} \&\>\> ${output_string}" | \
			xargs \
				-I input \
				-P $n_threads \
				bash -c "$gatk --java-options \"${java_options}\" input" \
				|| return 1
}

run_align_in_parallel() {
	local \
		aligner=$1 \
		parameter_file=$2 \
		align_options=$3 \
		n_threads=$4 \
		output_string=$5 
		sort_options=$6

	sed '/^#/d' ${parameter_file} \
		| parallel \
			--col-sep '\t' \
			echo -e "$align_options 2\>\> ${output_string} \| $samtools view -bS - 2\>\> ${output_string} \| $samtools sort $sort_options 2\>\> ${output_string}" \
			| xargs \
				-I input \
				-P $n_threads \
				bash -c "$aligner input" \
				|| return 1
}

set_up_tmps_and_logs() {
	local \
		cohort_log_file \
		id=$(random_id)
	
	if [[ -z ${inputs["log_prefix"]} ]]; then
		inputs["log_prefix"]=${log_dir}/${id}/ && mkdir ${inputs["log_prefix"]}
		# create cohort log file
		cohort_log_file=${inputs["log_prefix"]}${inputs["cohort_id"]}.log
		[[ -s $cohort_log_file ]] || > $cohort_log_file
		echo "==> output logs will be saved to ${inputs["log_prefix"]}"
	fi

	if [[ -z ${inputs["tmp_prefix"]} ]]; then
		inputs["tmp_prefix"]=${tmp_dir}/${id}_
		echo "==> temp. files will be saved to ${inputs["tmp_prefix"]}"
	fi
}

clean_temp() {
	[ ! -z ${inputs['tmp_prefix']} ] && rm ${inputs['tmp_prefix']}* 2> /dev/null

	return 0
}