#
#	NGSPIPELINE VCALL - pipeline module
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

#  call_variants()
#
#	* call variants using user-specified aligner
#	* aligner_id
#	* sw_id
#	* bam
#	* ref_id
#
call_variants() {

	dbsnpidx="dbsnp_${inputs['ref_id']}"

	if [ "${inputs['caller_id']}" = "gatk2_ug" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_prefix']}${inputs['cohort_id']}.${inputs['sample_id']}.log
		java -Xmx${inputs['maxmem']} -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} -jar ${gatk2}/GenomeAnalysisTK.jar -R ${b37_fasta} -T UnifiedGenotyper -glm BOTH --num_threads ${inputs['threads']} -I ${bam} --dbsnp ${!dbsnpidx} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "${inputs['caller_id']}" = "gatk3_hc" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		java -Xmx${inputs['maxmem']} -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} -jar ${gatk3}/GenomeAnalysisTK.jar -R ${b37_fasta} -T HaplotypeCaller -I ${bam} --dbsnp ${!dbsnpidx} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "${inputs['caller_id']}" = "gatk3_ug" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		java -Xmx${inputs['maxmem']} -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} -jar ${gatk3}/GenomeAnalysisTK.jar -R ${b37_fasta} -T UnifiedGenotyper -glm BOTH --num_threads ${inputs['threads']} -I ${bam} --dbsnp ${!dbsnpidx} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	
	elif [ "${inputs['caller_id']}" = "freebayes" ];then

		custom_call_2 _call_with_freebayes '  calling variants with freebayes...' || return 1

	elif [ "${inputs['caller_id']}" = "atlas_snp" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		ruby ${atlas_snp}/Atlas-SNP2.rb -i ${bam} -r ${b37_fasta} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf -n ${inputs['sample_id']} --Illumina
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "${inputs['caller_id']}" = "atlas_ind" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		ruby ${atlas_indel}/Atlas-Indel2.rb -b ${bam} -r ${b37_fasta} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf -I -s ${inputs['sample_id']}
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "${inputs['caller_id']}" = "glfsingle" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED piling up bam file for ${inputs['prefix']}" >> ${inputs['log_file']}
		${samtools_hybrid}/samtools-hybrid pileup -g -f ${b37_fasta} ${bam} > ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.glf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE piling up bam file for ${inputs['prefix']}" >> ${inputs['log_file']}
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		${glfsingle}/glfSingle -g ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.glf -b ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "${inputs['caller_id']}" = "ivc" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		cd ${work_dir}/${inputs['prefix']}
		cp ${ivc}/etc/ivc_config_default.ini ${work_dir}/${inputs['prefix']}/config.ini
		${ivc}/bin/configureWorkflow.pl --bam ${bam} --ref=${b37_fasta} --config=${work_dir}/${inputs['prefix']}/config.ini --output=${work_dir}/${inputs['prefix']}/${inputs['sample_id']}_${bam_type}
		cd ${work_dir}/${inputs['prefix']}/${inputs['sample_id']}_${bam_type} && make -j ${inputs['threads']}
		ln -s ${work_dir}/${inputs['prefix']}/${inputs['sample_id']}_${bam_type}/results/* ${work_dir}/${inputs['prefix']}
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	
	elif [ "${inputs['caller_id']}" = "samtools" ];then
		
		custom_call_2 _call_with_samtools '  calling variants with samtools...' || return 1

	elif [ "${inputs['caller_id']}" = "varscan" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		${samtools}/samtools mpileup -f ${b37_fasta} ${bam} | java -Xmx10g -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir} -jar ${varscan}/VarScan.v2.3.7.jar mpileup2cns --output-vcf 1 --variants > ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "${inputs['caller_id']}" = "platypus" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		python ${platypus}/Platypus.py callVariants --refFile ${b37_fasta} --bamFiles ${bam} --output ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_prefix']}${inputs['cohort_id']}.${inputs['sample_id']}.log
	fi

}

_call_with_freebayes() {
	local \
		input_bam_stage=$1 \
		log_file_string \
		line line_arr \
		sample_id \
		input_bam_files

	input_bam_files=()
	while read -r line; do
		[[ $line == "#"* ]] && continue
		line_arr=($line)
		sample_id=${line_arr[0]}

		input_bam_files+=(${bam_dir}/${inputs['prefix']}.${sample_id}.leftalignindels.bam)

	done < ${inputs['samples_list']}

	log_file_string="${inputs['log_prefix']}${inputs['cohort_id']}.log"

	echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${log_file_string}
	${freebayes} \
		-f ${inputs['ref']} \
		${input_bam_files[@]} \
		> ${vcf_dir}/${inputs['prefix']}.genotyped.g.vcf \
		2>> $log_file_string \
		|| return 1
	echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${log_file_string}
}

_call_with_samtools() {
	local \
		input_bam_stage=$1 \
		log_file_string \
		line line_arr \
		sample_id \
		input_bam_files
			
	input_bam_files=()
	while read -r line; do
		[[ $line == "#"* ]] && continue
		line_arr=($line)
		sample_id=${line_arr[0]}

		input_bam_files+=(${bam_dir}/${inputs['prefix']}.${sample_id}.leftalignindels.bam)

	done < ${inputs['samples_list']}

	log_file_string="${inputs['log_prefix']}${inputs['cohort_id']}.log"

	echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> $log_file_string
	${samtools} mpileup \
		-ugf ${inputs['ref']} \
		${input_bam_files[@]} \
		2>> $log_file_string \
		| ${bcftools} call \
			-vmO v \
			-o ${vcf_dir}/${inputs['prefix']}.genotyped.g.vcf \
			2>> $log_file_string \
			|| return 1
	echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> $log_file_string

}