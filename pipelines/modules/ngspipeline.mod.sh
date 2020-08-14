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

#  aligner()
#
# 	* sw_id=$1 - argument 1: the name of aligner
#	* ref_id=$2 - argument 2: the version of reference
#	* sample_id=$3 - argument 3: prefix for outputs
#	* fastq1=$4 - argument 4,5: input fastq files (gzipped)
#	* fastq2=$5
#
aligner() {

	echo "[${inputs['prefix']}][$(date "+%F %T")] +++ Start processing sample ${inputs['prefix']} with ${inputs['sw_id']}" >> ${inputs['log_file']}

	echo "[${inputs['prefix']}][$(date "+%F %T")] 1. STARTED mapping with ${inputs['prefix']}, convert, sort bam" >> ${inputs['log_file']}

	if [[ "${inputs['sw_id']}"="bwa" ]]; then
		${bwa} mem \
			-M \
			-t ${inputs['threads']} \
			-R "@RG\tID:${inputs['sample_id']}\tSM:${inputs['sample_id']}\tPL:Illumina" \
			${inputs['idx']} ${inputs['fastq1']} ${inputs['fastq2']} \
			| ${samtools} view \
				-bS - \
			| ${samtools} sort \
				-@ ${inputs['threads']} \
				-m ${inputs['maxmem']} \
				-o ${bam_dir}/${inputs['prefix']}.bam \
				-T ${tmp_dir}

	elif [ "${inputs['sw_id']}" = "gsnap" ];then
		${gsnap}/gsnap --gunzip -D ${!inputs['idx']} -d ${inputs['ref_id']} -t ${inputs['threads']} -A sam --npath=1 -B 5 --read-group-id ${inputs['sample_id']} --read-group-name ${inputs['sample_id']} --read-group-platform Illumina ${inputs['fastq1']} ${inputs['fastq2']} | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${inputs['threads']} -m ${inputs['maxmem']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	elif [ "${inputs['sw_id']}" = "bowtie2" ];then
		${bowtie2}/bowtie2 -p ${inputs['threads']} -x ${!inputs['idx']} -1 ${inputs['fastq1']} -2 ${inputs['fastq2']} --rg-id "${inputs['sample_id']}" --rg "SM:${inputs['sample_id']}" --rg "PL:Illumina" | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${inputs['threads']} -m ${inputs['maxmem']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	elif [ "${inputs['sw_id']}" = "soap2" ];then
		${soap2}/soap -p ${inputs['threads']} -a ${inputs['fastq1']} -b ${inputs['fastq2']} -D ${!inputs['idx']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_pe_output -2 ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_se_output -u ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_unmapped_output
		perl ${soap2sam} -p -s ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_se_output -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.sam -x ${inputs['sample_id']} ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_pe_output
		${samtools}/samtools view -bS ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.sam | ${samtools}/samtools sort -@ ${inputs['threads']} -m ${inputs['maxmem']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	elif [ "${inputs['sw_id']}" = "isaac" ];then
		isaac_fastq=${work_dir}/${inputs['prefix']}/isaac_fastq
		mkdir -p ${isaac_fastq}
		zcat ${inputs['fastq1']} | perl -pe 's/^\+ERR[0-9]+.[0-9]+ [0-9]+ length=101/\+/' | gzip - > ${isaac_fastq}/lane1_read1.fastq.gz
		zcat ${inputs['fastq2']} | perl -pe 's/^\+ERR[0-9]+.[0-9]+ [0-9]+ length=101/\+/' | gzip - > ${isaac_fastq}/lane1_read2.fastq.gz
		${isaac}/isaac-align --temp-directory ${temp_dir}/${inputs['prefix']} -r ${!inputs['idx']} -b ${isaac_fastq} -m 50 -j 6 --base-calls-format fastq-gz  --bam-header-tag "@RG\tID:${inputs['sample_id']}\tSM:${inputs['sample_id']}\tPL:Illumina" -o ${work_dir}/${inputs['prefix']} --verbosity 3 --input-parallel-load ${inputs['threads']} --temp-parallel-load ${inputs['threads']}
		ln -s ${work_dir}/${inputs['prefix']}/Projects/default/default/sorted.bam ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam
		ln -s ${work_dir}/${inputs['prefix']}/Projects/default/default/sorted.bam.bai ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam.bai
	
	elif [ "${inputs['sw_id']}" = "stampy" ];then
		python ${stampy}/stampy.py -g ${!inputs['idx']} -h ${!inputs['idx']} -t${inputs['threads']} --readgroup=ID:${inputs['sample_id']},SM:${inputs['sample_id']},PL:Illumina --bamkeepgoodreads -M ${bwa_bam} | ${samtools}/samtools view -bS - | ${samtools}/samtools sort -@ ${inputs['threads']} -m ${inputs['maxmem']} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	elif [ "${inputs['sw_id']}" = "novoalign" ];then
		${novoalign}/novoalign -c ${inputs['threads']} --mmapoff -t 20,3 --softclip 20 -d ${!inputs['idx']} -f ${inputs['fastq1']} ${inputs['fastq2']} -i 350 50 -k -o SAM "@RG\tID:${inputs['sample_id']}\tSM:${inputs['sample_id']}\tPL:Illumina" | ${samtools}/samtools view -Sb - | ${samtools}/samtools sort -@ ${inputs['threads']} -m 4G -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.bam -T ${tmp_dir}/${inputs['prefix']}
	
	fi

	echo "[${inputs['prefix']}][$(date "+%F %T")] 1. DONE mapping with ${inputs['prefix']}, convert, sort bam" >> ${inputs['log_file']}

	echo "[${inputs['prefix']}][$(date "+%F %T")] 2. STARTED mark duplicate/fixmate" >> ${inputs['log_file']}

	$gatk \
		--java-options "\
			-Xmx${inputs['maxmem']}	\
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']}" \
		MarkDuplicates \
			-I ${bam_dir}/${inputs['prefix']}.bam \
			-REMOVE_DUPLICATES true \
			-VALIDATION_STRINGENCY LENIENT \
			-AS true \
			-OUTPUT ${bam_dir}/${inputs['prefix']}.dedup.bam \
			-METRICS_FILE ${bam_dir}/${inputs['prefix']}.dedup.log \
			-CREATE_INDEX true

	$gatk \
		--java-options "\
			-Xmx${inputs['maxmem']}	\
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']}" \
		FixMateInformation \
			-I ${bam_dir}/${inputs['prefix']}.dedup.bam \
			-SORT_ORDER coordinate \
			-VALIDATION_STRINGENCY LENIENT \
			-CREATE_INDEX true \
			-O ${bam_dir}/${inputs['prefix']}.fixmate.bam

	echo "[${inputs['prefix']}][$(date "+%F %T")] 2. DONE mark duplicate/fixmate" >> ${inputs['log_file']}

	if [ "${recal_realign_on}" = "yes" ]||[ "${recal_realign_on}" = "both" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] 3. STARTED realign/recalibration" >> ${inputs['log_file']}

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3}/GenomeAnalysisTK.jar \
			-T RealignerTargetCreator \
			-nt ${inputs['threads']} \
			--known ${dbsnp_b37} \
			-R ${b37_fasta} \
			--filter_reads_with_N_cigar \
			-I ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_dedup_fixmate.bam \
			-o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.intervals

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3}/GenomeAnalysisTK.jar \
			-T IndelRealigner \
			-R ${b37_fasta} \
			--filter_reads_with_N_cigar \
			-I ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_dedup_fixmate.bam \
			-targetIntervals ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.intervals \
			--consensusDeterminationModel KNOWNS_ONLY \
			-known ${b37_1000g} \
			-LOD 0.4 \
			-o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_realign.bam

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3}/GenomeAnalysisTK.jar \
			-T BaseRecalibrator \
			--filter_reads_with_N_cigar \
			-I ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_realign.bam \
			-R ${b37_fasta} \
			-knownSites ${dbsnp_b37} \
			-o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.recal.table

		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3}/GenomeAnalysisTK.jar \
			-T PrintReads \
			--filter_reads_with_N_cigar \
			-I ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_realign.bam \
			-R ${b37_fasta} \
			-BQSR ${work_dir}/${inputs['prefix']}/${inputs['prefix']}.recal.table \
			-o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_recal_realign.bam

		echo "[${inputs['prefix']}][$(date "+%F %T")] 3. DONE realign/recalibration" >> ${inputs['log_file']}
	fi

	echo "[${inputs['prefix']}][$(date "+%F %T")] 4. STARTED left-align indels" >> ${inputs['log_file']}

	if [ "${recal_realign_on}" = "no" ]||[ "${recal_realign_on}" = "both" ];then
		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3}/GenomeAnalysisTK.jar \
			-R ${b37_fasta} \
			-T LeftAlignIndels \
			-I ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_dedup_fixmate.bam \
			-o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_org_leftalignindels.bam

	elif [ "${recal_realign_on}" = "yes" ]||[ "${recal_realign_on}" = "both" ];then
		java \
			-Xmx${inputs['maxmem']} \
			-XX:ParallelGCThreads=${inputs['threads']} \
			-Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} \
			-jar ${gatk3}/GenomeAnalysisTK.jar \
			-R ${b37_fasta} \
			-T LeftAlignIndels \
			-I ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_recal_realign.bam \
			-o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_recal_leftalignindels.bam
	fi

	echo "[${preifx}][$(date "+%F %T")] 4. DONE left-align indels" >> ${inputs['log_file']}
	echo "[${inputs['prefix']}][$(date "+%F %T")] +++ End processing sample ${inputs['prefix']} with ${inputs['sw_id']}" >> ${inputs['log_file']}

}

#  vcaller()
#
#	* call variants...
#
vcaller() {
	local \
		aligner_id=$1 \
		sw_id=$2 \
		bam=$3 \
		ref_id=$4

	if [ $# -gt 4 ];then
		local bam_type=$5
	else
		local bam_type="recal"
	fi

	local prefix=${inputs['sample_id']}_${aligner_id}_${inputs['sw_id']}_${inputs['ref_id']}

	if [ ! -d  ${tmp_dir}/${inputs['prefix']} ];then
		mkdir ${tmp_dir}/${inputs['prefix']}
	fi

	if [ ! -d ${work_dir}/${inputs['prefix']} ];then
		mkdir ${work_dir}/${inputs['prefix']}
		mkdir ${work_dir}/${inputs['prefix']}/log
	fi

	dbsnpidx="dbsnp_${inputs['ref_id']}"
	LOG_FILE="${work_dir}/${inputs['prefix']}/log/${inputs['prefix']}.log"

	echo "[${inputs['prefix']}][$(date "+%F %T")] +++ Start processing sample ${inputs['sample_id']} with ${inputs['sw_id']}" >> ${inputs['log_file']}

	if [ "$sw_id" = "gatk2_ug" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		java -Xmx${inputs['maxmem']} -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} -jar ${gatk2}/GenomeAnalysisTK.jar -R ${b37_fasta} -T UnifiedGenotyper -glm BOTH --num_threads ${inputs['threads']} -I ${bam} --dbsnp ${!dbsnpidx} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 200
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "gatk3_hc" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		java -Xmx${inputs['maxmem']} -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} -jar ${gatk3}/GenomeAnalysisTK.jar -R ${b37_fasta} -T HaplotypeCaller -I ${bam} --dbsnp ${!dbsnpidx} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "gatk3_ug" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		java -Xmx${inputs['maxmem']} -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir}/${inputs['prefix']} -jar ${gatk3}/GenomeAnalysisTK.jar -R ${b37_fasta} -T UnifiedGenotyper -glm BOTH --num_threads ${inputs['threads']} -I ${bam} --dbsnp ${!dbsnpidx} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "freebayes" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		${freebayes}/freebayes -f ${b37_fasta} ${bam} > ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "atlas_snp" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		ruby ${atlas_snp}/Atlas-SNP2.rb -i ${bam} -r ${b37_fasta} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf -n ${inputs['sample_id']} --Illumina
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "atlas_ind" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		ruby ${atlas_indel}/Atlas-Indel2.rb -b ${bam} -r ${b37_fasta} -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf -I -s ${inputs['sample_id']}
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "glfsingle" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED piling up bam file for ${inputs['prefix']}" >> ${inputs['log_file']}
		${samtools_hybrid}/samtools-hybrid pileup -g -f ${b37_fasta} ${bam} > ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.glf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE piling up bam file for ${inputs['prefix']}" >> ${inputs['log_file']}
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		${glfsingle}/glfSingle -g ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.glf -b ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "ivc" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		cd ${work_dir}/${inputs['prefix']}
		cp ${ivc}/etc/ivc_config_default.ini ${work_dir}/${inputs['prefix']}/config.ini
		${ivc}/bin/configureWorkflow.pl --bam ${bam} --ref=${b37_fasta} --config=${work_dir}/${inputs['prefix']}/config.ini --output=${work_dir}/${inputs['prefix']}/${inputs['sample_id']}_${bam_type}
		cd ${work_dir}/${inputs['prefix']}/${inputs['sample_id']}_${bam_type} && make -j ${inputs['threads']}
		ln -s ${work_dir}/${inputs['prefix']}/${inputs['sample_id']}_${bam_type}/results/* ${work_dir}/${inputs['prefix']}
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "samtools" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		${samtools}/samtools mpileup -ugf ${b37_fasta} ${bam} | ${bcftools}/bcftools call -vmO v -o ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "varscan" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		${samtools}/samtools mpileup -f ${b37_fasta} ${bam} | java -Xmx10g -XX:ParallelGCThreads=${inputs['threads']} -Djava.io.tmpdir=${tmp_dir} -jar ${varscan}/VarScan.v2.3.7.jar mpileup2cns --output-vcf 1 --variants > ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	elif [ "$sw_id" = "platypus" ];then
		echo "[${inputs['prefix']}][$(date "+%F %T")] STARTED calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
		python ${platypus}/Platypus.py callVariants --refFile ${b37_fasta} --bamFiles ${bam} --output ${work_dir}/${inputs['prefix']}/${inputs['prefix']}_${bam_type}.vcf
		echo "[${inputs['prefix']}][$(date "+%F %T")] DONE calling variants with ${inputs['prefix']}" >> ${inputs['log_file']}
	fi

	echo "[${inputs['prefix']}][$(date "+%F %T")] +++ End processing sample ${inputs['sample_id']} with ${inputs['sw_id']}" >> ${inputs['log_file']}
}