#
#  SETUP - setup script 
#		for sequence-analysis-pipelines library
#
#	* run this script first before the sap.sh wrapper,
#	* but after some tools have been installed. 
#	*
#	* You may re-run this script any time another
#	* tool has been installed to update your instance
#	* of the sap library.
#	*
#	* basic usage:
#	* ./setup.sh <study-type> <reference-id>
#
#	* example:
#	* ./setup.sh wgs lambda-virus
#
#	Jack Morrice
#
##########################################################
#!/bin/bash

source includes/utilities.sh

workflow() {
	local argv=("$@")

	declare -A inputs=( ['log_file']=$(random_id)_setup.log ['study_type']=${argv[0]} ['ref_label']=${argv[1]})
	declare -A direcs=( ['pro_dir']=$(pwd) )
	declare -A tools

	echo "==> any output logs will be saved to ${inputs['log_file']}"

	custom_call initialize_direcs_hash "setting the project sub-directory locations hash..."

	custom_call initialize_tools_hash "setting the tool path hash..."

	# format data directory
	custom_call format_data_directory "formatting data directory for analysis..."

	# check which tools have been installed
	custom_call check_installed_tools "checking which tools have been installed..."

	# create or update includes/locations.sh
	custom_call create_locations_file " writing includes/locations.sh"

	# build all installed tool indices
	custom_call build_tool_indices "building indices for all installed tools..."
}


#
#  tasks
#
initialize_direcs_hash() {
	direcs['src_dir']=${direcs['pro_dir']}/source
	direcs['tls_dir']=${direcs['pro_dir']}/tools
	direcs['dat_dir']=${direcs['pro_dir']}/data
	direcs['sam_dir']=${direcs['dat_dir']}/sam
	direcs['bam_dir']=${direcs['dat_dir']}/bam
	direcs['vcf_dir']=${direcs['dat_dir']}/vcf
	direcs['rds_dir']=${direcs['dat_dir']}/reads
	direcs['ref_dir']=${direcs['dat_dir']}/references
	direcs['fqc_dir']=${direcs['dat_dir']}/fastqc
	direcs['tmp_dir']=${direcs['dat_dir']}/temp
	direcs['log_dir']=${direcs['dat_dir']}/logs
	direcs['med_dir']=${direcs['dat_dir']}/media
	direcs['pip_dir']=${direcs['pro_dir']}/pipelines
	direcs['pin_dir']=${direcs['pro_dir']}/inputs
	direcs['tst_dir']=${direcs['pro_dir']}/test-suites
}

initialize_tools_hash() {

	# set defaults
	tools['bowtie2']=NULL
	tools['bowtie2_build']=NULL
	tools['samtools']=NULL
	tools['bcftools']=NULL
	tools['bgzip']=NULL
	tools['bwa']=NULL
	tools['gmap_build']=NULL
	tools['gsnap']=NULL
	tools['stampy']=NULL
	tools['gatk']=NULL
	tools['fastqc']=NULL
	tools['trimmomatic']=NULL
	tools['freebayes']=NULL

	## these tools are inlcuded in the github repo (no need for installation)
	tools['simulate']=${direcs['tls_dir']}/simulate-0.1/simulate.pl
	tools['jaccard']=${direcs['tls_dir']}/jaccard/jaccard.py
	tools['generatejson']=${direcs['tls_dir']}/jaccard/generate-json.py
	tools['benchmarkR']=${direcs['tls_dir']}/jaccard/benchmark.R

	# update values from tool list json:
	_update_tool_value_from_list 'bowtie2'
	_update_tool_value_from_list 'bowtie2_build'
	_update_tool_value_from_list 'samtools'
	_update_tool_value_from_list 'bcftools'
	_update_tool_value_from_list 'bgzip'
	_update_tool_value_from_list 'bwa'
	_update_tool_value_from_list 'gmap_build'
	_update_tool_value_from_list 'gsnap'
	_update_tool_value_from_list 'stampy'
	_update_tool_value_from_list 'gatk'
	_update_tool_value_from_list 'fastqc'
	_update_tool_value_from_list 'trimmomatic'
	_update_tool_value_from_list 'freebayes'
}

_update_tool_value_from_list() {
	local key=$1
	value_from_json ${direcs['tls_dir']}/tool-list.json '.'${key} tools[${key}]
	if [[ ! ${tools[${key}]} == "NULL" ]]; then
		tools[${key}]=${direcs['tls_dir']}/${tools[${key}]}
	fi
}

format_data_directory() {
	local status=0

	[[ ! -d ${direcs['bam_dir']} ]] &&  { mkdir ${direcs['bam_dir']}; echo "  created ${direcs['bam_dir']}"; }
	[[ ! -d ${direcs['fqc_dir']} ]] &&  { mkdir ${direcs['fqc_dir']}; echo "  created ${direcs['fqc_dir']}"; }
	[[ ! -d ${direcs['log_dir']} ]] &&  { mkdir ${direcs['log_dir']}; echo "  created ${direcs['log_dir']}"; }
	[[ ! -d ${direcs['med_dir']} ]] &&  { mkdir ${direcs['med_dir']}; echo "  created ${direcs['med_dir']}"; }
	[[ ! -d ${direcs['rds_dir']} ]] &&  { mkdir ${direcs['rds_dir']}; echo "  created ${direcs['rds_dir']}"; }
	[[ ! -d ${direcs['sam_dir']} ]] &&  { mkdir ${direcs['sam_dir']}; echo "  created ${direcs['sam_dir']}"; }
	[[ ! -d ${direcs['tmp_dir']} ]] &&  { mkdir ${direcs['tmp_dir']}; echo "  created ${direcs['tmp_dir']}"; }
	[[ ! -d ${direcs['vcf_dir']} ]] &&  { mkdir ${direcs['vcf_dir']}; echo "  created ${direcs['vcf_dir']}"; }

	return $status
}

check_installed_tools() {
	local status=0

	if [[ ! ${tools['bowtie2']} == "NULL" ]]; then
		printf '  checking bowtie2...'
		if _check_bowtie2 ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['samtools']} == "NULL" ]]; then
		printf '  checking samtools...'
		if _check_samtools ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['bcftools']} == "NULL" ]]; then
		printf '  checking bcftools...'
		if _check_bcftools ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['bgzip']} == "NULL" ]]; then
		printf '  checking bgzip...'
		if _check_bgzip ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['bwa']} == "NULL" ]]; then
		printf '  checking bwa...'
		if _check_bwa ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['gsnap']} == "NULL" ]]; then
		printf '  checking gsnap...'
		if _check_gsnap ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['stampy']} == "NULL" ]]; then
		printf '  checking stampy...'
		if _check_stampy ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['gatk']} == "NULL" ]]; then
		printf '  checking gatk...'
		if _check_gatk ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['fastqc']} == "NULL" ]]; then
		printf '  checking fastqc...'
		if _check_fastqc ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['trimmomatic']} == "NULL" ]]; then
		printf '  checking trimmomatic...'
		if _check_trimmomatic ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	if [[ ! ${tools['freebayes']} == "NULL" ]]; then
		printf '  checking freebayes...'
		if _check_freebayes ; then
			echo "...installed"
		else
			echo "...not installed!"
			status=1
		fi
	fi

	return $status
}


_check_bowtie2() { ${tools['bowtie2']} --help &> /dev/null || return 1; }
_check_samtools() {	${tools['samtools']} --help &> /dev/null || return 1; }
_check_bcftools() {	${tools['bcftools']} --help &> /dev/null || return 1; }
_check_bgzip() { ${tools['bgzip']} --help &> /dev/null || return 1; }
_check_bwa() { man ${tools['bwa']}.1 &> /dev/null || return 1; }
_check_gsnap() { ${tools['gsnap']} --help &> /dev/null || return 1; }
_check_stampy() { python ${tools['stampy']} --help &> /dev/null || return 1; }
_check_gatk() { ${tools['gatk']} --help &> /dev/null || return 1; }
_check_fastqc() { ${tools['fastqc']} --help &> /dev/null || return 1; }
_check_trimmomatic() { java -jar ${tools['trimmomatic']} -version &> /dev/null || return 1; }
_check_freebayes() { ${tools['freebayes']} --help &> /dev/null || return 1; }



build_tool_indices() {
	local status=0

	if [[ ! ${tools['bwa']} == "NULL" ]]; then
		printf '  building bwa index...'
		if _build_bwa_index ; then
			echo "...done"
		else
			echo "...failed!"
			status=1
		fi
	fi

	if [[ ! ${tools['gsnap']} == "NULL" ]]; then
		printf '  building gmap index...'
		if _build_gmap_index ; then
			echo "...done"
		else
			echo "...failed!"
			status=1
		fi
	fi

	if [[ ! ${tools['bowtie2']} == "NULL" ]]; then
		printf '  building bowtie2 index...'
		if _build_bowtie2_index ; then
			echo "...done"
		else
			echo "...failed!"
			status=1
		fi
	fi

	if [[ ! ${tools['stampy']} == "NULL" ]]; then
		printf '  building stampy index and hash...'
		if _build_stampy_index ; then
			echo "...done"
		else
			echo "...failed!"
			status=1
		fi
	fi

	if [[ ! ${tools['samtools']} == "NULL" ]]; then
		printf '  building samtools index...'
		if _build_samtools_index ; then
			echo "...done"
		else
			echo "...failed!"
			status=1
		fi
	fi

	if [[ ! ${tools['gatk']} == "NULL" ]]; then
		printf '  building gatk dict...'
		if _build_gatk_dict ; then
			echo "...done"
		else
			echo "...failed!"
			status=1
		fi
	fi

	return $status
}

_build_bwa_index() {

	# option -a specifies the algorithm, with arguments:
	# 
	# 1. 'is' is for short sequences (like virus genomes)
	# 2. 'bwtsw' is for long sequences (like human genomes)

	if [[ ! -f "${direcs['ref_dir']}/${inputs['ref_label']}/bwa.${inputs['ref_label']}.amb" ]]; then
		${tools['bwa']} index \
			-p ${direcs['ref_dir']}/${inputs['ref_label']}/bwa.${inputs['ref_label']} \
			-a is \
			${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}.fa.gz \
			&>> ${inputs['log_file']}
	else
		printf "skipping as index exists"
	fi
}

_build_gmap_index() {
	# this index is used by both gmap and gsnap
	if [[ ! -d "${direcs['ref_dir']}/${inputs['ref_label']}/gmap.${inputs['ref_label']}" ]]; then
		${tools['gmap_build']} \
			-D ${direcs['ref_dir']}/${inputs['ref_label']} \
			-d "gmap.${inputs['ref_label']}" \
			-g \
			${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}'.fa.gz' \
			&>> ${inputs['log_file']}
	else
		printf "skipping as index exists"
	fi
}

_build_bowtie2_index() {
	if [[ ! -f "${direcs['ref_dir']}/${inputs['ref_label']}/bowtie2.${inputs['ref_label']}.1.bt2l" ]]; then
		${tools['bowtie2_build']} \
			--large-index \
			${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}'.fa.gz' \
			${direcs['ref_dir']}/${inputs['ref_label']}/bowtie2.${inputs['ref_label']} \
			&>> ${inputs['log_file']}
	else
		printf "skipping as index exists"
	fi
}

_build_stampy_index() {

	# build index
	if [[ ! -f ${direcs['ref_dir']}/${inputs['ref_label']}/stampy.${inputs['ref_label']}.stidx ]]; then
		python ${tools['stampy']} \
			-G ${direcs['ref_dir']}/${inputs['ref_label']}/stampy.${inputs['ref_label']} \
			${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}'.fa.gz' \
			&>> ${inputs['log_file']}
	else
		printf "skipping as index "
	fi

	# build hash
	if [[ ! -f ${direcs['ref_dir']}/${inputs['ref_label']}/stampy.${inputs['ref_label']}.sthash ]]; then
		python ${tools['stampy']} \
			-g ${direcs['ref_dir']}/${inputs['ref_label']}/stampy.${inputs['ref_label']} \
			-H ${direcs['ref_dir']}/${inputs['ref_label']}/stampy.${inputs['ref_label']} \
			&>> ${inputs['log_file']}
	else
		printf "and hash exists"
	fi

}

_build_samtools_index() {

	if [[ ! -f "${direcs['ref_dir']}/${inputs['ref_label']}/samtools.${inputs['ref_label']}.fai" ]]; then
		${tools['samtools']} faidx \
			${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}'.fa.gz' \
			> ${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}.fai \
			2>> ${inputs['log_file']}
	else
		printf "skipping as index exists"
	fi

}

_build_gatk_dict() {

	if [[ ! -f "${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}.dict" ]]; then
		${tools['gatk']} CreateSequenceDictionary \
			-R ${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}'.fa.gz' \
			-O ${direcs['ref_dir']}/${inputs['ref_label']}/${inputs['ref_label']}'.dict' \
			&>> ${inputs['log_file']}
	else
		printf "skipping as dictionary exists"
	fi
}


create_locations_file() {
	local x

	echo '''#
#  LOCATIONS - bash config file
#
#	* file containing all path definitions
#	* for directories and executable tools
#	* used by the variant-caller package.
#	* 
#	* Use the <sap-path>/setup.sh script to 
#	* modify this file. Any direct changes
#	* to the file itself will be overwritten
#	* by another call of setup.sh.
#
#	Jack Morrice
#
################################################
	''' > includes/locations.sh

	echo '# project path variables' >> includes/locations.sh
	for x in "${!direcs[@]}"
	do 
		printf "%s=%s\n" "$x" "${direcs[$x]}" >> includes/locations.sh
	done
	echo >> includes/locations.sh

	echo '# tool executable paths' >> includes/locations.sh
	for x in "${!tools[@]}"
	do 
		if [[ ! "${tools[$x]}" == "NULL" ]]; then
			printf "%s=%s\n" "$x" "${tools[$x]}" >> includes/locations.sh
		fi
	done
}

#
#  run workflow
#
workflow "$@"