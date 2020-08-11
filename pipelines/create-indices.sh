#
#  CREATE INDICES - setup script for pipelines
#
#	* given an input reference sequence, build
#	* all necessary indexes for the different 
#	* aliner and variant caller algorithms to
#	* be used by the pipelines in this library.
#
#	Jack Morrice
#
##################################################
#!/bin/bash

source ../includes/locations.sh
source ../includes/utilities.sh

source ${pip_dir}/modules/checkparams.mod.sh


workflow() { 
	local argv=("$@")

	declare -A inputs=( ["input_json"]=${argv[0]} )

	custom_call check_input_json "checking a pipeline input json file was provided..."

	custom_call initialize_inputs_hash "initializing pipeline input parameter values..."

	custom_call build_bwa_index "building bwa alignment index..."

	custom_call build_reference_index "building reference index..."

	custom_call build_reference_dict "building reference dictionary"
}

#
#  tasks
#
initialize_inputs_hash() {
	local status=0

	# 1. set default parameter values
	printf '  setting default parameter values...'
	inputs["ref"]=NULL
	echo '...done'

	# 2. update values from input json file
	printf '  updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.ref_base'	   inputs["ref_base"] && inputs["ref"]=${ref_dir}/${inputs["ref_base"]}
	echo '...done'

	# 3. check input parameters make sense
	printf '  checking that parameter values make sense...'
	check_ref || status=1
	[[ $status == 0 ]] && echo '...done'

	return $status
}


build_bwa_index() {
	local \
		ref_file \
		idx_file

	# option -a specifies the algorithm, with arguments:
	# 
	# 1. 'is' is for short sequences (like virus genomes)
	# 2. 'bwtsw' is for long sequences (like human genomes)

	ref_file=${inputs['ref']}
	idx_file=${ref_file/.fa/.bwa}

	if [[ ! -f "${idx_file}.amb" ]]; then
		$bwa index \
			-p ${ref_file/.fa/.bwa} \
			-a is \
			$ref_file
	else
		echo "skipping bwa alignment index build as index already exists"
	fi
	
}

build_reference_index() {
	local \
		ref_file \
		idx_file

	# the gatk haplotype variant caller requires a dictionary and an index 
	# are both built from the reference. this function builds the index 
	# which I think is different from the alignment index?
	#
	# http://www.htslib.org/doc/samtools-faidx.html
	# https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference

	ref_file=${inputs['ref']}
	idx_file=${ref_file/.fa/.fai}

	if [[ ! -f "$idx_file" ]]; then
		$samtools faidx \
			$ref_file \
			> $idx_file
	else
		echo "skipping reference index build as this index already exists"
	fi

}

build_reference_dict() {
	local \
		ref_file \
		dict_file

	# the gatk haplotype variant caller requires a dictionary and an index 
	# are both built from the reference. this function builds the dictionary 
	# which I think is different from the alignment index?
	#
	# https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
	# https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-

	ref_file=${inputs['ref']}
	dict_file=${ref_file/.fa/dict}

	if [[ ! -f "$dict_file" ]]; then
		$gatk CreateSequenceDictionary \
			-R $ref_file \
			-O $dict_file
	else
		echo "skipping reference dictionary build as this dictionary already exists"
	fi

}


#
#  run workflow
#
workflow "$@"

