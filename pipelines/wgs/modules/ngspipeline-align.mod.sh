#
#	NGSPIPELINE ALIGN - pipeline module
#
#	* This module contains all subroutines
#	* related to alignment for the ngspipeline.
#	*
#	* ngspipeline was adapted from a paper in Nature
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

#  align_reads_to_reference()
#
#	* this function is the only one that
#	* should be called outside this module. 
#	* the other functions are all subroutines
#	* of this one. 
#
align_reads_to_reference() {

	custom_call_2 _align_cohort_reads "  aligning with ${inputs['aligner_id']}..." || return 1

	custom_call_2 _mark_cohort_duplicates '  marking duplicates...' || return 1

	custom_call_2 _fix_cohort_mate_info '  fixing mate information...' || return 1

	if [ "${inputs['recal_realign_on']}" = "yes" ]; then
		custom_call_2 _perform_cohort_realignments '  performing realignments...' || return 1
	fi

	if [ "${inputs['recal_realign_on']}" = "no" ]; then
		custom_call_2 _left_align_cohort_indels '  left-aligning indels...' fixmate || return 1
	elif [ "${inputs['recal_realign_on']}" = "yes" ]; then
		custom_call_2 _left_align_cohort_indels '  left-aligning indels...' realign || return 1
	fi
}

#  subroutines for align_reads_to_reference()
#
#	* these functions are all subroutines
#	* and not to be used outside this module. 
#
_align_cohort_reads() {
	cohort_call _align_sample_reads "mapping reads with ${inputs['aligner_id']}" || return 1
}

_mark_cohort_duplicates() {
	local \
		java_options \
		option_string \
		log_file_string

	java_options="-Xmx${inputs['maxmem']}	\
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']}"

	option_string="\
		MarkDuplicates \
			-I ${bam_dir}/${inputs['prefix']}.{1}.bam \
			-REMOVE_DUPLICATES true \
			-VALIDATION_STRINGENCY LENIENT \
			-AS true \
			-OUTPUT ${bam_dir}/${inputs['prefix']}.{1}.dedup.bam \
			-METRICS_FILE ${bam_dir}/${inputs['prefix']}.{1}.dedup.log \
			-CREATE_INDEX true"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{1}.log"

	run_gatk_in_parallel \
		"${java_options}" \
		${inputs['samples_list']} \
		"${option_string}" \
		${inputs['threads']} \
		${log_file_string}

}

_fix_cohort_mate_info() {
	cohort_call _fix_sample_mate_info "fixing mate information in sample sequence" || return 1
}

_left_align_cohort_indels() {
	local input_bam_stage=$1
	cohort_call _left_align_sample_indels "left-aligning indels in sample sequence" $input_bam_stage || return 1
}

_perform_cohort_realignments() {
	cohort_call _perform_sample_realignments "realigning sample sequences with BaseRecalibrator etc" || return 1
}

cohort_call() {
	local \
		routine=$1 \
		log_msg=$2 \
		log_file_string \
		line line_arr \
		sample_id \
		fastq1 fastq2
		shift; shift

	while read -r line; do
		[[ $line == "#"* ]] && continue
		line_arr=($line)
		sample_id=${line_arr[0]}
		fastq1=${rds_dir}/${line_arr[1]}
		fastq2=${rds_dir}/${line_arr[2]}

		log_file_string="${inputs['log_prefix']}${inputs['cohort_id']}.${sample_id}.log"

		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED ${log_msg}" >> $log_file_string

		$routine $log_file_string $sample_id $fastq1 $fastq2 "$@" || return 1

		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE ${log_msg}" >> ${inputs['log_prefix']}${inputs['cohort_id']}.${sample_id}.log

	done < ${inputs['samples_list']}
}

_align_sample_reads() {
	local \
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

	if [ ${inputs['aligner_id']} = "bwa" ]; then

		_align_with_bwa $log_file_string $sample_id $fastq1 $fastq2 || return 1

	elif [ ${inputs['aligner_id']} = "gsnap" ]; then

		_align_with_gsnap $log_file_string $sample_id $fastq1 $fastq2 || return 1
	
	elif [ ${inputs['aligner_id']} = "bowtie2" ]; then

		_align_with_bowtie2 $log_file_string $sample_id $fastq1 $fastq2 || return 1
	
	elif [ "${inputs['aligner_id']}" = "soap2" ];then
		${soap2}/soap -p ${inputs['threads']} -a ${inputs['fastq1']} -b ${inputs['fastq2']} -D ${!inputs['idx']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_pe_output -2 ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_se_output -u ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_unmapped_output
		perl ${soap2sam} -p -s ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_se_output -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.sam -x ${inputs['sample_id']} ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_pe_output
		${samtools}/samtools view -bS ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.sam | ${samtools}/samtools sort -@ ${inputs['threads']} -m ${inputs['maxmem']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	elif [ "${inputs['aligner_id']}" = "isaac" ];then
		isaac_fastq=${work_dir}/${inputs['prefix']}/isaac_fastq
		mkdir -p ${isaac_fastq}
		zcat ${inputs['fastq1']} | perl -pe 's/^\+ERR[0-9]+.[0-9]+ [0-9]+ length=101/\+/' | gzip - > ${isaac_fastq}/lane1_read1.fastq.gz
		zcat ${inputs['fastq2']} | perl -pe 's/^\+ERR[0-9]+.[0-9]+ [0-9]+ length=101/\+/' | gzip - > ${isaac_fastq}/lane1_read2.fastq.gz
		${isaac}/isaac-align --temp-directory ${temp_dir}/${inputs['prefix']} -r ${!inputs['idx']} -b ${isaac_fastq} -m 50 -j 6 --base-calls-format fastq-gz  --bam-header-tag "@RG\tID:${inputs['sample_id']}\tSM:${inputs['sample_id']}\tPL:Illumina" -o ${work_dir}/${inputs['prefix']} --verbosity 3 --input-parallel-load ${inputs['threads']} --temp-parallel-load ${inputs['threads']}
		ln -s ${work_dir}/${inputs['prefix']}/Projects/default/default/sorted.bam ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam
		ln -s ${work_dir}/${inputs['prefix']}/Projects/default/default/sorted.bam.bai ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam.bai
	
	elif [ "${inputs['aligner_id']}" = "stampy" ]; then

		_align_with_stampy $log_file_string $sample_id $fastq1 $fastq2 || return 1
	
	# novoalign is not free, so this option is depricated:
	elif [ "${inputs['aligner_id']}" = "novoalign" ];then
		${novoalign}/novoalign -c ${inputs['threads']} --mmapoff -t 20,3 --softclip 20 -d ${!inputs['idx']} -f ${inputs['fastq1']} ${inputs['fastq2']} -i 350 50 -k -o SAM "@RG\tID:${inputs['sample_id']}\tSM:${inputs['sample_id']}\tPL:Illumina" | ${samtools}/samtools view -Sb - | ${samtools}/samtools sort -@ ${inputs['threads']} -m 4G -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	fi
}

_align_with_bwa() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

	${bwa} mem \
		-M \
		-t ${inputs['threads']} \
		-R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:Illumina" \
		${ref_dir}/bwa.${inputs['ref_base']%.fa} ${fastq1} ${fastq2} \
		2>> $log_file_string \
		| ${samtools} view \
			-bS - \
			2>> $log_file_string \
			| ${samtools} sort \
				-@ ${inputs['threads']} \
				-m ${inputs['maxmem']} \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.bam \
				-T ${tmp_dir} \
				2>> $log_file_string
}

_align_with_gsnap() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

	${gsnap} \
		--gunzip \
		--dir ${ref_dir} \
		--db gmap.${inputs['ref_base']%.fa} \
		-t ${inputs['threads']} \
		-A sam \
		--npath=1 \
		-B 5 \
		--read-group-id ${inputs['sample_id']} \
		--read-group-name ${inputs['sample_id']} \
		--read-group-platform Illumina \
		${inputs['fastq1']} ${inputs['fastq2']} \
		2>> $log_file_string \
		| ${samtools} view \
			-bS - \
			2>> $log_file_string \
			| ${samtools} sort \
				-@ ${inputs['threads']} \
				-m ${inputs['maxmem']} \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.bam \
				-T ${tmp_dir} \
				2>> $log_file_string

}

_align_with_bowtie2() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

	${bowtie2} \
		-p ${inputs['threads']} \
		-x ${ref_dir}/bowtie2.${inputs['ref_base']%.fa} \
		-1 ${inputs['fastq1']} \
		-2 ${inputs['fastq2']} \
		--rg-id "${inputs['sample_id']}" \
		--rg "SM:${inputs['sample_id']}" \
		--rg "PL:Illumina" \
		2>> $log_file_string \
		| ${samtools} view \
			-bS - \
			2>> $log_file_string \
			| ${samtools} sort \
				-@ ${inputs['threads']} \
				-m ${inputs['maxmem']} \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.bam \
				-T ${tmp_dir} \
				2>> $log_file_string
}

_align_with_stampy() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

	${bwa} mem \
		-M \
		-t ${inputs['threads']} \
		-R "@RG\tID:${inputs['sample_id']}\tSM:${inputs['sample_id']}\tPL:Illumina" \
		${ref_dir}/bwa.${inputs['ref_base']%.fa} ${inputs['fastq1']} ${inputs['fastq2']} \
		2>> $log_file_string \
		| ${samtools} view \
			-bS - \
			2>> $log_file_string \
			| ${samtools} sort \
				-@ ${inputs['threads']} \
				-m ${inputs['maxmem']} \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.bwa.bam \
				-T ${tmp_dir} \
				2>> $log_file_string

	python ${stampy} \
		-g ${ref_dir}/stampy.${inputs['ref_base']%.fa} \
		-h ${ref_dir}/stampy.${inputs['ref_base']%.fa} \
		-t${inputs['threads']} \
		--readgroup=ID:${inputs['sample_id']},SM:${inputs['sample_id']},PL:Illumina \
		--bamkeepgoodreads -M ${bam_dir}/${inputs['prefix']}.${sample_id}.bwa.bam \
		2>> $log_file_string \
		| ${samtools} view \
			-bS - \
			2>> $log_file_string \
			| ${samtools} sort \
				-@ ${inputs['threads']} \
				-m ${inputs['maxmem']} \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.bam \
				-T ${tmp_dir} \
				2>> $log_file_string
}

_mark_sample_duplicates() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

	$gatk \
		--java-options "\
			-Xmx${inputs['maxmem']}	\
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']}" \
		MarkDuplicates \
			-I ${bam_dir}/${inputs['prefix']}.${sample_id}.bam \
			-REMOVE_DUPLICATES true \
			-VALIDATION_STRINGENCY LENIENT \
			-AS true \
			-OUTPUT ${bam_dir}/${inputs['prefix']}.${sample_id}.dedup.bam \
			-METRICS_FILE ${bam_dir}/${inputs['prefix']}.${sample_id}.dedup.log \
			-CREATE_INDEX true \
		&>> $log_file_string

}

_fix_sample_mate_info() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4
	
	$gatk \
		--java-options "\
			-Xmx${inputs['maxmem']}	\
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']}" \
		FixMateInformation \
			-I ${bam_dir}/${inputs['prefix']}.${sample_id}.dedup.bam \
			-SORT_ORDER coordinate \
			-VALIDATION_STRINGENCY LENIENT \
			-CREATE_INDEX true \
			-O ${bam_dir}/${inputs['prefix']}.${sample_id}.fixmate.bam \
		&>> $log_file_string
}

_left_align_sample_indels() {
	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4 \
		input_bam_stage=$5

	$gatk \
		--java-options "\
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']}" \
		LeftAlignIndels \
			-R ${inputs['ref']} \
			-I ${bam_dir}/${inputs['prefix']}.${sample_id}.${input_bam_stage}.bam \
			-O ${bam_dir}/${inputs['prefix']}.${sample_id}.leftalignindels.bam \
		&>> $log_file_string
}

_perform_sample_realignments() {

	# this function has not yet been tested, 
	# and currently *does not work*.

	local 
		log_file_string=$1 \
		sample_id=$2 \
		fastq1=$3 \
		fastq2=$4

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3} \
				-T RealignerTargetCreator \
				-nt ${inputs['threads']} \
				--known ${dbsnp_b37} \
				-R ${inputs['ref']} \
				--filter_reads_with_N_cigar \
				-I ${bam_dir}/${inputs['prefix']}.${sample_id}.fixmate.bam \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.intervals \
			2>> $log_file_string

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3} \
				-T IndelRealigner \
				-R ${inputs['ref']} \
				--filter_reads_with_N_cigar \
				-I ${bam_dir}/${inputs['prefix']}.${sample_id}.fixmate.bam \
				-targetIntervals ${bam_dir}/${inputs['prefix']}.intervals \
				--consensusDeterminationModel KNOWNS_ONLY \
				#-known ${b37_1000g} \
				-LOD 0.4 \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.realign.bam
			2>> $log_file_string

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3} \
				-T BaseRecalibrator \
				--filter_reads_with_N_cigar \
				-I ${bam_dir}/${inputs['prefix']}.${sample_id}.realign.bam \
				-R ${inputs['ref']} \
				#-knownSites ${dbsnp_b37} \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.recal.table
			2>> $log_file_string

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3} \
				-T PrintReads \
				--filter_reads_with_N_cigar \
				-I ${bam_dir}/${inputs['prefix']}.${sample_id}.realign.bam \
				-R ${inputs['ref']} \
				-BQSR ${bam_dir}/${inputs['prefix']}.${sample_id}.recal.table \
				-o ${bam_dir}/${inputs['prefix']}.${sample_id}.recal.realign.bam
			2>> $log_file_string
}
