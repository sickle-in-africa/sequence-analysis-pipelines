#
#  function library for geneMapNGS pipeline
#
#	Kevin Esoh
#
#	contents:
#		1. command functions
#			*these functions are called at the top level
#			by the pipeline scripts*
#		2. check functions
#			*these functions check that the necessary 
#			files exist for the command functions to run
#			properly (called by the command functions 
#			themselves; not called at pipeline level).*
#		3. prep functions
#			*these functions prepare the user's project
#			directory for the command functions to run
#			(again, called by the command functions 
#			themselves; not called at pipeline level)*
#
#############################################################

function start() {
    echo -e """
               ===================================================================
               \e[38;5;43mGeneMAP NGS Pipeline			          GeneMAP (c) 2020\e[0m
               -------------------------------------------------------------------
               Argument:            Parameter
               --------             --------
               Fastq/SAM/BAM path:  $dname
               reference:           $ref
	       	   dbsnp:		        $dbsnp
               leading:             $leadx
               trailing:            $trailx
               PED file:            $ped
               threads:             $t
               adapter:             $adap
               sample file:         $meta
               BAM list:            $blist
               GVCF list:           $glist
               known sites:         $ks
               outFile:             ${out}.vcf.gz
               ===================================================================
               Starting NGS Pipeline. Please wait...
    """
}

initialize_inputs_hash() {
	local status=0

	# 1. set default parameter values
	printf '  setting default parameter values...'
	inputs["meta"]=NULL
	inputs["cohort_id"]=NULL
	inputs["threads"]=1
	inputs["trim_minlen"]=36
	inputs["trim_adap"]=NULL
	inputs["trim_leadx"]=0
	inputs["trim_trailx"]=0
	inputs["trim_av_qual_min"]=0
	inputs["ref"]=NULL
	inputs["blist"]=NULL
	inputs["glist"]=NULL
	inputs["ped"]=NULL
	inputs["dbsnp"]=NULL
	echo '...done'

	# 2. update parameters with arguments from the input json file
	printf '  updating with arguments from input json file...'
	value_from_json ${inputs["input_json"]} '.cohort_id'   inputs["cohort_id"]
	value_from_json ${inputs["input_json"]} '.meta_base'   inputs["meta_base"]
	if [[ ! "${inputs["cohort_id"]}" == NULL ]]; then
		inputs["meta"]=${rds_dir}/${inputs["cohort_id"]}.txt
	elif [[ ! "${inputs["meta_base"]}" == NULL ]]; then
		inputs["meta"]=${rds_dir}/${inputs["meta_base"]}
	fi
	value_from_json ${inputs["input_json"]} '.threads'     inputs["threads"]
	value_from_json ${inputs["input_json"]} '.trim_minlen' inputs["trim_minlen"]
	value_from_json ${inputs["input_json"]} '.trim_leadx'  inputs["trim_leadx"]
	value_from_json ${inputs["input_json"]} '.trim_trailx' inputs["trim_trailx"]
	value_from_json ${inputs["input_json"]} '.trim_av_qual_min' inputs["trim_av_qual_min"]
	value_from_json ${inputs["input_json"]} '.ref_base'	   inputs["ref_base"] && inputs["ref"]=${ref_dir}/${inputs["ref_base"]}
	value_from_json ${inputs["input_json"]} '.blist' 	   inputs["blist"]
	value_from_json ${inputs["input_json"]} '.glist' 	   inputs["glist"]
	value_from_json ${inputs["input_json"]} '.ped' 	   	   inputs["ped"]
	value_from_json ${inputs["input_json"]} '.dbsnp' 	   inputs["dbsnp"]
	echo '...done'

	# 3. check that inputs make sense
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
	## still need to write checks for the other input parameters...
	[[ $status == 0 ]] && echo '...done'

	# 4. set up output logging and temporary files
	set_up_tmps_and_logs || { echo 'seting up temp and log directory failed'; status=1; }

	return $status
}

#
#  1. command functions
#
#  summary:
#	+ fq(), pfq()
#	+ trim(), ptrim()
#	+ bmap(), pbmap()
#	+ gmap(), pgmap()
#	+ indelreal(), pindelreal()
#	+ bqsr(), pbqsr()
#	+ emit_gvcfs(), pemit_gvcfs()
#	+ combinegvcfs()
#	+ genogvcfs()
#	+ varcall()
#	+ vqsr()
#	+ bcfcall()
#	

#  fq()
#
#	*calls FastQC on read files in series, parallel*
#
#	input:
#		+ $meta
#	tools required:
#		+ $fastqc
#		+ gnu parallel (for pfq() only)
#	outputs:
#		+ fastq quality report files in ${fqc_dir}
#		+ forward_reverse.txt in ${tmp_dir} 
#
function fq() {
	local \
		option_string \
		log_file_string

	option_string="\
		-t ${inputs["threads"]} \
		-o ${fqc_dir}/ \
		${rds_dir}/{1} ${rds_dir}/{2}"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  checking fastq quality with FASTQC..."
	run_in_parallel \
		$fastqc \
		${inputs["meta"]} \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

#  tim()
#
#	*clips unwanted reads in fastq files with trimmomatic*
#	inputs:
#		+ $meta
#		+ forward_reverse.txt in ${tmp_dir} 
#	tools required:
#		+ $trimmomatic
#	outputs:
#		+ trimmed (and zipped) read (fastq) files
#
function trim() {
	local \
		option_string \
		log_file_string

	option_string="-jar ${trimmomatic} PE \
		-threads ${inputs["threads"]} \
		-basein ${rds_dir}/${inputs["cohort_id"]}.{3}.raw_1P.fq \
		-baseout ${rds_dir}/${inputs["cohort_id"]}.{3}.trimmed.fq \
		LEADING:${inputs["trim_leadx"]} \
		TRAILING:${inputs["trim_trailx"]} \
		SLIDINGWINDOW:5:${inputs["trim_av_qual_min"]} \
		MINLEN:${inputs["trim_minlen"]}"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  trimming reads with trimmomatic..."
	run_in_parallel \
		"java" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }	
}

#  bmap(), pbmap()
#
#	* Mapping/Alignment (BWA)
#	* this function does three things in sequence:
#	*	1. aligns all reads to the reference sequence with bwa mem
#	*		 to create sam files...
#	*	2. ...converts the sam files to bam files with $samtools view...
#	*	3. ...and then sorts the bam files into 'mapped' bam files,
#	*		and removes the unsorted (unmapped bam files).
#
#	inputs:
#		+ $t
#		+ $ref
#
#	tools required:
#		+ $bwa
#		+ $samtools
#
#	outputs:
#		+ aligned and sorted (mapped) sequences for each sample, 
#			saved as compressed bam files
#
function bmap() {
	
	prep_reads_list 'trimmed' || return 1

	_aligning_reads_to_reference || return 1

	_converting_sam_to_bam || return 1

	_sorting_bam_files || return 1

	_add_read_group_info || return 1

	_mark_duplicates || return 1

	_validate_bam_files || return 1
}

#  gmap()
#
#	* GATKv4 BWA Mapping/Alignment
#	* its odd: only the parallel version marks duplicates
#	* and indexes the final bam files. Is this a mistake?
# 
#	inputs:
#		+ $meta
#		+ $t
#		+ $ref
#		+ a bwa index for $ref
#		+ a gatk dictionary for $ref
#
#	tools required:
#
#	outputs:
#		+ [gmap()]  ${bam_dir}/<sample_id>.merged.bam (mapped bam files for each sample)
#
function gmap() {

	prep_reads_list 'trimmed' || return 1

	_converting_fastq_to_sam || return 1

	_aligning_reads_to_reference || return 1

	_converting_sam_to_bam || return 1

	_sorting_bam_files || return 1

	_merge_bam_alignments || return 1

	_mark_duplicates || return 1

	_validate_bam_files || return 1

}

_converting_fastq_to_sam() {
	local \
		option_string \
		log_file_string \
		status=0

	option_string="FastqToSam \
		-F1 ${rds_dir}/{2} \
		-F2 ${rds_dir}/{3} \
		-SM {1} \
		-PL {4} \
		-RG {1} \
		-O ${bam_dir}/${inputs["cohort_id"]}.{1}.unmapped.bam"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{1}.log"

	printf "  converting fastq files to unmapped bam..."
	run_in_parallel \
		"$gatk" \
		"${inputs['tmp_prefix']}reads.list" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_aligning_reads_to_reference() {
	local \
		option_string \
		log_file_string \
		status=0

	option_string="mem \
		-t ${inputs['threads']} \
		${inputs['ref']} ${rds_dir}/{2} ${rds_dir}/{3} \
		-o ${sam_dir}/${inputs["cohort_id"]}.{1}.sam"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{1}.log"

	printf "  aligning reads to reference with bwa..."
	run_in_parallel \
		"$bwa" \
		"${inputs['tmp_prefix']}reads.list" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_converting_sam_to_bam() {
	local \
		option_string \
		status=0

	option_string="view \
		-h ${sam_dir}/${inputs["cohort_id"]}.{3}.sam \
		-O BAM \
		-o ${bam_dir}/${inputs["cohort_id"]}.{3}.bam"

	log_file_sting="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  converting mapped sam files to bam..."
	run_in_parallel \
		"$samtools" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_sting}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_sorting_bam_files() {
	local \
		option_string \
		status=0

	option_string="sort \
		-O BAM \
		--reference ${inputs["ref"]} \
		-@ ${inputs["threads"]} \
		-o ${bam_dir}/${inputs["cohort_id"]}.{3}.sorted.bam \
		${bam_dir}/${inputs["cohort_id"]}.{3}.bam"

	log_file_sting="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  sorting bam files..."
	run_in_parallel \
		"$samtools" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_sting}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_merge_bam_alignments() {
	local \
		option_string \
		status=0

	option_string="MergeBamAlignment \
		-R ${inputs["ref"]} \
		-UNMAPPED ${bam_dir}/${inputs["cohort_id"]}.{3}.unmapped.bam \
		-ALIGNED ${bam_dir}/${inputs["cohort_id"]}.{3}.sorted.bam \
		-O ${bam_dir}/${inputs["cohort_id"]}.{3}.merged.bam"

	log_file_sting="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  merging aligned and unaligned bam files..."
	run_in_parallel \
		"$gatk" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_sting}" \
		"${log_file_sting}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_add_read_group_info() {
	# this function may be useful when we don't use
	# the merge bam alignmnet step
	local \
		option_string \
		status=0

	option_string="AddOrReplaceReadGroups \
		-I ${bam_dir}/${inputs["cohort_id"]}.{3}.sorted.bam \
		-O ${bam_dir}/${inputs["cohort_id"]}.{3}.merged.bam \
		-SM {3} \
		-PL {4} \
		-PU {3} \
		-LB {3}"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  adding read group information to bam files..."
	run_in_parallel \
		"$gatk" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status

}

_mark_duplicates() {
	local \
		input_bam_stage='merged' \
		sample_id \
		line \
		option_string \
		status=0

	option_string="MarkDuplicates \
		-I ${bam_dir}/${inputs["cohort_id"]}.{3}.${input_bam_stage}.bam \
		-O ${bam_dir}/${inputs["cohort_id"]}.{3}.marked.bam \
		-M ${bam_dir}/${inputs["cohort_id"]}.{3}.marked_dup_metrics.txt"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  marking duplicates..."
	run_in_parallel \
		"$gatk" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_validate_bam_files() {
	local \
		option_string \
		status=0

	option_string="ValidateSamFile \
		-I ${bam_dir}/${inputs["cohort_id"]}.{3}.marked.bam \
		-R ${inputs["ref"]} \
		--TMP_DIR ${tmp_dir}/ \
		-M SUMMARY"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  validating bam files..."
	run_in_parallel \
		"$gatk" \
		"${inputs["meta"]}" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }

	return $status
}


#  indelreal(), pindelreal()
#
#	* Indel Realignment 
#	* (This will not be run if GATK v4 is used since HaplotypeCaller 
#	* essentially does local rearrangements).
#
function indelreal() {
	local \

	prep_bam_list ${inputs['tmp_prefix']} 'marked' || return 1

	_create_realigner_targets ${inputs['tmp_prefix']}bam.list || return 1

	_realign_indels || return 1

	prep_bam_list ${inputs['tmp_prefix']} 'realigned' || return 1

	_index_bam_files ${inputs['tmp_prefix']}bam.list || return 1

	return 0

       id="bam.list"
       n=$((50/$t))
       mkdir -p realigned
       while read -r line; do
            $gatk \
            	-T RealignerTargetCreator \
            	-R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "rtc.ks.txt" -a -s "rtc.ks.txt" ]; then cat rtc.ks.txt; fi; fi) \
            	-I aligned/${line} \
            	-o realigned/${line/.bam/.intervals}
            $gatk \
            	-T IndelRealigner \
            	-R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "ir.ks.txt" -a -s "ir.ks.txt" ]; then cat ir.ks.txt; fi; fi) \
            	-I aligned/${line} \
            	-targetIntervals realigned/${line/.bam/.intervals} \
            	-o realigned/${line/.bam/.realigned.bam}
            $samtools index \
            	-b aligned/${line/.bam/.realigned.bam}
       done < ${id}
       if [ -e "ir.ks.txt" -o -e "rtc.ks.txt" ]; then rm ir.ks.txt || rm rtc.ks.txt; fi
}

_create_realigner_targets() {
	local \
		bamlist=$1 \
		option_string \
		status=0

	java -jar $gatk3 -T RealignerTargetCreator
	exit 0

	option_string="-jar ${gatk3} \
		-T RealignerTargetCreator \
		-R ${inputs["ref"]} \
		-I ${bam_dir}/{2} \
		-o ${bam_dir}/{2}.intervals"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{1}.log"

	printf "  creating realigner targets..."
	run_in_parallel \
		"java" \
		"$bamlist" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_realign_indels() {
	local \
		bamlist=$1 \
		option_string \
		status=0

	option_string="-jar $gatk3 \
		-T IndelRealigner \
		-R ${inputs["ref"]} \
		-I ${bam_dir}/{2} \
		-targetIntervals ${bam_dir}/{2}.intervals \
		-o ${bam_dir}/{3}.realigned.bam"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{1}.log"

	printf "  realigning indels..."
	run_in_parallel \
		"java" \
		"$bamlist" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; status=1; } \
		&& { echo "...done."; }

	return $status
}




function pindelreal() {
       check_ref; check_bamlist
       #--- Indel realignment (This requires GATKv3.x. Point to your installation of it in 'gatk3_esoh' above)
       ###  According to GATK Best Practices, this step is not necessary in the new pipeline, as HaplotypeCaller does a good job  ###
       #awk '{print $1,$2,$3,$4}' ${meta} > metadat.txt
       id="bam.list"
       n=$((50/$t))
       mkdir -p realigned
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo -T RealignerTargetCreator -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "rtc.ks.txt" -a -s "rtc.ks.txt" ]; then cat rtc.ks.txt; fi; fi) -I aligned/{1}.bam -o realigned/{1}.intervals | xargs -I input -P$n sh -c "$gatk3 input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo -T IndelRealigner -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "ir.ks.txt" -a -s "ir.ks.txt" ]; then cat ir.ks.txt; fi; fi) -I aligned/{1}_mkdups.bam -targetIntervals realigned/{1}.intervals -o realigned/{1}.realigned.bam | xargs -I input -P$n sh -c "$gatk3 input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo index -b aligned/{1}.realigned.bam | xargs -I input -P$n sh -c "samtools input"
       if [ -e "ir.ks.txt" -o -e "rtc.ks.txt" ]; then rm ir.ks.txt || rm rtc.ks.txt; fi
}

#  bqsr(), pbqsr()
#
#	* Base Quality Score Recallibration (BQSR)
#
#	inputs:
#		+ bam list (optional)
#		+ The input read data whose base quality scores need to be assessed.
#		+ A database of known polymorphic sites to skip over.
#	tools:
#		+ gatk 4.0 (BaseRecalibrator, Apply BQSR)
#	outputs:
#		+ A GATK Report file with many tables
function bqsr() {

	check_ref || return 1
	check_gatk_dict || return 1
	check_bamlist $tmp_prefix || return 1

	id=${tmp_prefix}bam.list
	n=$((50/$t))

	while read -r line; do
		echo $line
		return 1
		$gatk BaseRecalibrator \
			-I ${dname}/${line/.bam/} \
			-R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "bqsr.ks.txt" -a -s "bqsr.ks.txt" ]; then rm rtc.ks.txt ir.ks.txt; cat bqsr.ks.txt; fi; fi) \
			-O bqsr/${line/.bam/_recal_data.table} \
			|| { echo "BaseRecalibrator step failed"; return 1; }
		$gatk ApplyBQSR \
			-R $ref 
			-I ${dname}/${line/.bam/} 
			--bqsr-recal-file bqsr/${line/.bam/_recal_data.table} \
			-O aligned/${line/.bam/.mapped.bam} \
			|| { echo "ApplyBQSR step failed"; return 1; }
	done < ${id}
	if [ -e "${id}" ]; then rm ${id}; fi

	# remove all temporary files
	[ ! -z "${tmp_prefix}" ] && rm ${tmp_prefix}*

}

function pbqsr() {
       if [[ "$dname" == NULL ]]; then
          echo -e "\n\e[38;5;1mERROR\e[0m: -p,--path not provided! Please specify path to BAM files"; 1>&2;
          return 1;
       fi
       check_ref; check_gatk_dict; check_bamlist
       id="$blist"
       n=$((50/$t))
       mkdir -p bqsr; mkdir -p aligned
       #--- Base quality score recalibration (BQSR)
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo BaseRecalibrator -I ${dname}/{1/}.bam -R $ref $(if [[ $ks != NULL ]]; then check_sites; if [ -e "bqsr.ks.txt" -a -s "bqsr.ks.txt" ]; then rm rtc.ks.txt ir.ks.txt; cat bqsr.ks.txt; fi; fi) -O bqsr/{1/}_recal_data.table | xargs -I input -P$n sh -c "gatk input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' echo ApplyBQSR -R $ref -I ${dname}/{1/}.bam --bqsr-recal-file bqsr/{1/}_recal_data.table -O aligned/{1/}.mapped.bam | xargs -I input -P$n sh -c "gatk input"
       sed 's/.bam//g' ${id} | parallel --col-sep ' ' rm ${dname}/{1/}.bam
       if [ -e "${id}" ]; then rm ${id}; fi
}

#  varcall
#
#	* Variant Calling with GATK (Single Cohort Joint) in Serial
#	* This function takes analysis ready bam files, and:
#	* 1. indexes them
#	* 2. calls variants with the gatk haplotype caller
#	* 3. ...
#	* the bam files are the output of the bqsr step of the 
#	* pipeline. 
#	inputs:
#
#	tools required:
#
#	outputs:
#
function varcall() {

	prep_bam_list 'marked' || return 1

	_index_bam_files || return 1

	_call_sample_variants || return 1

	prep_gvcf_list 'raw' || return 1

	_combine_sample_gvcfs  || return 1

	_genotype_combined_gvcf || return 1

}

_index_bam_files() {
	local \
		option_string \
		status=0

	option_string="index \
		-b \
		-@ ${inputs["threads"]} \
		${bam_dir}/{2}"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  indexing bam files..."
	run_in_parallel \
		"$samtools" \
		"${inputs['tmp_prefix']}bam.list" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; }

	return $status
}

_call_sample_variants() {
	local \
		option_string \
		log_file_string

	option_string="HaplotypeCaller \
		-R ${inputs["ref"]} \
		-I ${bam_dir}/{2} \
		-O ${vcf_dir}/${inputs['cohort_id']}.{1}.raw.g.vcf \
		$(if [[ ${inputs["ped"]} != NULL ]]; then echo -ped ${inputs["ped"]}; fi) \
		$(if [[ ${inputs["dbsnp"]} != NULL ]]; then echo --dbsnp ${inputs["dbsnp"]}; fi) \
		--lenient true \
		-ERC GVCF"

	log_file_string="${inputs["log_prefix"]}${inputs["cohort_id"]}.{3}.log"

	printf "  calling sample variants with gatk HaplotypeCaller..."
	run_in_parallel \
		"$gatk" \
		"${inputs['tmp_prefix']}bam.list" \
		"${option_string}" \
		"${inputs["threads"]}" \
		"${log_file_string}" \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

_combine_sample_gvcfs() {

	printf "  combining sample gvcf files..."
	$gatk CombineGVCFs \
		-R ${inputs["ref"]} \
		--arguments_file ${inputs['tmp_prefix']}gvcf.list \
		$(if [[ ${inputs["dbsnp"]} != NULL ]]; then echo --dbsnp ${inputs["dbsnp"]}; fi) \
		$(if [[ ${inputs["ped"]} != NULL ]]; then echo -ped {inputs["ped"]}; fi) \
		-O ${vcf_dir}/${inputs["cohort_id"]}.combined.g.vcf \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

_genotype_combined_gvcf() {

	printf "  genotyping cohort's combined gvcf file..."
	$gatk GenotypeGVCFs \
		-R ${inputs["ref"]} \
		-V ${vcf_dir}/${inputs["cohort_id"]}.combined.g.vcf \
		$(if [[ ${inputs["dbsnp"]} != NULL ]]; then echo --dbsnp ${inputs["dbsnp"]}; fi) \
		$(if [[ ${inputs["ped"]} != NULL ]]; then echo -ped {inputs["ped"]}; fi) \
		-O ${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }

}


function vqsr() {
   vcf=$1; op=$2
 
   #SNPs
   gatk VariantRecalibrator \
      -R /mnt/lustre/groups/CBBI1243/KEVIN/db/ucsc.hg19.fasta \
      -V ${vcf} \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/hapmap_3.3.hg19.sites.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/1000G_phase1.indels.hg19.sites.vcf.gz \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/1000G_omni2.5.hg19.sites.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /mnt/lustre/groups/CBBI1243/KEVIN/db/dbsnp_138.hg19.vcf.gz \
      -an QD \
      -an MQ \
      -an MQRankSum \
      -an ReadPosRankSum \
      -an FS \
      -an SOR \
      -an InbreedingCoeff \
      -mode BOTH \
      -O ${op}.recal \
      --tranches-file ${op}.tranches \
      --rscript-file ${op}.plots.R
   
   #Apply
   gatk ApplyVQSR \
       -V ${vcf} \
       --recal-file ${op}.recal \
       -O ${op}.vqsr-filtered.vcf.gz 

}

#  bcfcall
#
#	* Variant Calling With bcftools
#	inputs:
#		+ list of analysis-ready (marked) bam files
#	tools:
#		+ bcftools (mpileup, index, call, view)
#	outputs:
#		+ a gvcf file containing variants jointly called across the cohort
#
function bcfcall() {

	prep_bam_list 'marked' || return 1

	cat ${inputs['tmp_prefix']}bam.list \
		| awk -v d=${bam_dir}/ '{print d$2}' \
		> ${inputs['tmp_prefix']}bam.inputs

	_pile_up_sample_gvcfs || return 1

	_index_gvcf_file 'piledup' || return 1

	_joint_call_variants_bcftools || return 1

	_index_gvcf_file 'genotyped' || return 1

	_unzip_gvcf_file 'genotyped' || return 1

}

_pile_up_sample_gvcfs() {

	printf "  piling up sample gvcf files before joint calling..."
	$bcftools mpileup \
		--min-MQ 1 \
		--thread ${inputs["threads"]} \
		-f ${inputs["ref"]} \
		-Oz \
		-o ${vcf_dir}/${inputs["cohort_id"]}.piledup.g.vcf.gz \
		-b ${inputs['tmp_prefix']}bam.inputs \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

_index_gvcf_file() {
	local input_gvcf_stage=$1

	printf "  indexing cohort piledup gvcf file..."
	$bcftools index -f -t ${vcf_dir}/${inputs["cohort_id"]}.${input_gvcf_stage}.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

_joint_call_variants_bcftools() {

	printf "  joint-calling variants for cohort with bcftools..."
	$bcftools call \
		-mv \
		--threads ${inputs["threads"]} \
		-Oz \
		-o ${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf.gz \
		${vcf_dir}/${inputs["cohort_id"]}.piledup.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

_unzip_gvcf_file() {

	printf "  unzipping the joint-called cohort gvcf file..."
	$bcftools view \
		-o ${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf \
		${vcf_dir}/${inputs["cohort_id"]}.genotyped.g.vcf.gz \
		&>> ${inputs["log_prefix"]}${inputs["cohort_id"]}.log \
		|| { echo "...failed!"; return 1; } \
		&& { echo "...done."; return 0; }
}

#
#  2. check functions
#
#	summary:
#	+ check_sites()
#	+ check_adapter()
#	+ check_sample()
#	+ check_fq()
#	+ check_bamlist()
#	+ check_gvcflist()
#

#--- Check References and their indexes


#--- Check additional [optional] references (Known sites) for IndelRealignment, BQSR, and VQSR
function check_sites() {
    nks=$(echo $ks | sed 's/,/ --known /g')
    echo "--known $nks" > rtc.ks.txt
    nks=$(echo $ks | sed 's/,/ -known /g')
    echo "-known $nks" > ir.ks.txt
    nks=$(echo $ks | sed 's/,/ --known-sites /g')
    echo "--known-sites $nks" > bqsr.ks.txt
}

#--- Check trimmomatic adapters
function check_adapter() {
    function warning() {
        echo -e """\e[38;5;3mWARNING\e[0m: The adapter was not found! Make sure it is present in the current directory""" 1>&2;
        echo -e """\e[38;5;6m===>\e[0m Attempting to trim without adapter. Press \e[38;5;6mCTRL+C\e[0m to stop\n""" 1>&2;
        return 1;

    }
    case "$(echo "$adap" | tr [:lower:] [:upper:])" in
        NP) if [[ -e "NexteraPE-PE.fa" ]]; then echo ILLUMINACLIP:NexteraPE-PE.fa:2:30:10; else warning; fi ;;
        T3U) if [[ -e "TruSeq3-PE-2.fa" ]]; then echo ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10;  else warning; fi ;;
        T2P) if [[ -e "TruSeq2-PE.fa" ]]; then echo ILLUMINACLIP:TruSeq2-PE.fa:2:30:10;  else warning; fi ;;
        T3P) if [[ -e "TruSeq3-PE.fa" ]]; then echo ILLUMINACLIP:TruSeq3-PE.fa:2:30:10;  else warning; fi ;;
        T2S) if [[ -e "TruSeq2-SE.fa" ]]; then echo ILLUMINACLIP:TruSeq2-SE.fa:2:30:10;  else warning; fi ;;
        T3S) if [[ -e "TruSeq3-SE.fa" ]]; then echo ILLUMINACLIP:TruSeq3-SE.fa:2:30:10;  else warning; fi ;;
         *) echo -e """\e[38;5;3mWARNING\e[0m: No such adapter '$adap'! Type --help for usage\n\e[38;5;6m===>\e[0m Attempting to trim without adapter. Press \e[38;5;6mCTRL+C\e[0m to stop\n""" 1>&2; return 1; ;;
    esac
}

#--- Check sample file
function check_sample() {
	local \
		line \
		line_array \
		sample_id \
		sample_log_file

	# check meta variabe was set
	if [[ "${inputs["meta"]}" == NULL ]]; then
		echo -e "\e[38;5;1mERROR\e[0m: -s,--sample_list not provided! Exiting..."; 1>&2;
		return 1
	elif [ -f ${inputs["meta"]} -a -s ${inputs["meta"]} ]; then
		for i in $(awk '{print $1}' ${inputs["meta"]}); do
			[ ! $i == "#"* ] && continue
			if [ ! -f ${rds_dir}/$i ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' was not found in the directory '${rds_dir}'.\nPlease specify the path with -p or --path or check that the files in the path are the same in the sample list" 1>&2;
				return 1;
			elif [ -f ${rds_dir}/${i} -a ! -s ${rds_dir}/${i} ]; then
				echo -e "\e[38;5;1mERROR\e[0m: '${i}' may be empty. Please check and correct '${rds_dir}'." 1>&2;
				return 1;
		   fi
		done
	elif [ -f ${inputs["meta"]} -a ! -s ${inputs["meta"]} ]; then
		echo -e "\e[38;5;1mERROR\e[0m: '$meta' seems to be empty! Please check and correct." 1>&2;
	fi

	# check read files in $meta file exist, and create log files if they do not exist
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			sample_id=${line_array[2]}
			[[ -s ${rds_dir}/"${line_array[0]}" ]] || { echo "Error: ${line_array[0]} does not exist or is empty"; return 1; }
			[[ -s ${rds_dir}/"${line_array[1]}" ]] || { echo "Error: ${line_array[0]} does not exist or is empty"; return 1; }


		fi
	done < ${inputs["meta"]}

}

#--- Prepare input for BQSR
function check_gvcflist() {
if [[ ( "$glist" == NULL ) ]]; then
   if [ -f "gvcf.list" ]; then
      rm gvcf.list;
   fi;
   for i in *.gvcf*; do
      if [[ ( -f ${i} ) && ( -s ${i} ) ]]; then  # if gvcf files exist in the current directory and are not empty
         basename -a $(ls $i) >> gvcf.list;
      elif [[ -d vcall ]]; then # if a directory exists called vcall
         for j in vcall/*.gvcf*; do
             if [[ ( -f ${j} ) && ( -s ${j} ) ]]; then # if gvcf files exist in the vcall directory and are not empty
                basename -a $(ls $j) >> gvcf.list;
             fi;
         done;
      else
         echo -e "\n\e[38;5;1mERROR\e[0m: Please check that there are GVCF files in the path $dname\n" 1>&2;
         #return 1;
      fi;
   done;
   if [ -f "gvcf.list" -a -s "gvcf.list" ]; then
      echo -e "\n\e[38;5;6mNOTE\e[0m: $(cat gvcf.list | wc -l) GVCF file(s) counted in '$(readlink -f $(dirname $(cat gvcf.list | head -1)))/' and will be used! Press CTRL+C to stop\n";
      sleep 1;
   fi;
fi
}


#  3. prep functions
#
#	summary:
#	+ prep_fq()
#	+ prep_trim()
#	+ prep_map()
#	+ prep_bam_list()
#	+ prep_gvcf_list()
#

#  prep_fq
#
#	* create temporary files used by fq()
#
function prep_fq() {
	local \
		line

	#--- Make input files from forward/reverse runs or SAM/BAM files
	> ${inputs['tmp_prefix']}fwd.txt
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			echo -e "${line_array[3]}\t${line_array[0]}" >> ${inputs['tmp_prefix']}fwd.txt
		fi
	done < ${inputs["meta"]}

	> ${inputs['tmp_prefix']}rev.txt
	while IFS= read -r line; do
		if [[ ! $line == "#"* ]]; then
			line_array=( $line )
			echo -e "${line_array[3]}\t${line_array[1]}" >> ${inputs['tmp_prefix']}rev.txt
		fi
	done < ${inputs["meta"]}

	if [[ ! -s "${inputs['tmp_prefix']}rev.txt" ]]; then
		cp ${inputs['tmp_prefix']}fwd.txt ${inputs['tmp_prefix']}forward_reverse.txt
		awk -v d="${rds_dir}/" '{print $1 d$2}' ${inputs['tmp_prefix']}forward_reverse.txt > ${inputs['tmp_prefix']}fastq.input.txt
	else
		paste ${inputs['tmp_prefix']}fwd.txt ${inputs['tmp_prefix']}rev.txt | awk '{print $1,$2}' > ${inputs['tmp_prefix']}forward_reverse.txt
		awk -v d="${rds_dir}/" '{print $1 d$2,d$3}' ${inputs['tmp_prefix']}forward_reverse.txt > ${inputs['tmp_prefix']}fastq.input.txt
	fi
	rm ${inputs['tmp_prefix']}fwd.txt ${inputs['tmp_prefix']}rev.txt

}

function prep_trim() {
    
    prep_fq ${inputs['tmp_prefix']} || return 1

    #--- On checking for fastq files above, we checked for SAM/BAM as well. If the function picked SAM/BAM, we definitely wanna spill errors since we can't trim SAM/BAM here
    for i in $(awk '{print $1}' ${inputs['tmp_prefix']}forward_reverse.txt | head -1); do
        if [[ ( ${i} == *.sam ) || ( ${i} == *.sam.gz ) || ( "${i}" == *.bam ) ]]; then
           echo -e "\n\e[38;5;1mERROR\e[0m: No fastq/SAM/BAM file found in the specified location: '${rds_dir}'\nPlease specify path to Fastq/SAM/BAM files using -p or --path\n" 1>&2;
           return 1;
           rm ${inputs['tmp_prefix']}forward_reverse.txt ${inputs['tmp_prefix']}fastq.input.txt
        fi
    done
    awk -v d="${rds_dir}/" '{print d$1,d$2, d$1".paired_fp.fq.gz", d$1".unpaired_fu.fq.gz", d$2".paired_rp.fq.gz", d$2".unpaired_ru.fq.gz"}' ${inputs['tmp_prefix']}forward_reverse.txt > ${inputs['tmp_prefix']}trim.input.txt
	#rm forward_reverse.txt fastq.input.txt;
}

#--- Prepare alignment/mapping input
function prep_map() {
	tmp_prefix=$1
	check_ref || return 1
	check_fq $tmp_prefix || return 1
	prep_trim $tmp_prefix || return 1

	if [ -e "${tmp_prefix}forward_reverse.txt" -a -s "${tmp_prefix}forward_reverse.txt" ]; then
		awk -v i="${rds_dir}/" -v o="${sam_dir}/" '{print i$1,i$2,"-o",o$1".sam"}' ${tmp_prefix}forward_reverse.txt > ${tmp_prefix}align.input.txt
	else
		echo -e "\n\e[38;5;1mERROR\e[0m: Please check that there are fastq files in the path...\n"
	fi
}

function prep_reads_list() {
	local \
		input_reads_stage=$1

	> ${inputs['tmp_prefix']}reads.list
	while read -r line; do
		[[ "$line" == "#"* ]] && continue
		line_array=( $line )
		sample_id=${line_array[2]}
		reads_prefix="${inputs[cohort_id]}.${sample_id}.${input_reads_stage}"
		platform=${line_array[3]}

		echo -e "${sample_id}\t${reads_prefix}_1P.fq\t${reads_prefix}_2P.fq\t${platform}" >> ${inputs['tmp_prefix']}reads.list
	done < ${inputs["meta"]}
}

#  prep_bam_list
#
#  Prepare input for BQSR, varcall
#	* check that a list of bam files
#	* has been provided by the user
#	* via $blist, 
#	* or create one from the sample list
#	* via $meta
# 
function prep_bam_list() {
	#  * if the user has specified an input bam list 
	#  * then use this one, otherwise, create one
	#  * using the sample IDs in the $inputs["meta"] file.
	local \
		input_bam_stage=$1 \
		samples \
		i

	if [[ ! ${inputs["blist"]}==NULL ]]; then
		cat ${inputs["blist"]} > ${inputs['tmp_prefix']}bam.list 
		return 0
	else
		samples=$(sed '/^#/d' ${inputs["meta"]} | awk '{print $3}')
		for i in ${samples[@]}; do
			echo -e "$i\t${inputs["cohort_id"]}.$i.${input_bam_stage}.bam"
		done > ${inputs['tmp_prefix']}bam.list
	fi

	[[ -s ${inputs['tmp_prefix']}bam.list ]] || { echo "\n\e[38;5;1mERROR\e[0m: no bam files were found, please check and try again"; return 1; }

}

function prep_gvcf_list() {
	#  * if the user has specified an input gvcf list 
	#  * then use this one, otherwise, create one
	#  * using the sample IDs in the $inputs["meta"] file.
	local \
		input_gvcf_stage=$1 \
		samples \
		i

	if [[ ! ${inputs["glist"]}==NULL ]]; then
		echo "using user's input gvcf list $blist"
		cat ${inputs["glist"]} > ${inputs['tmp_prefix']}gvcf.list 
		return 0
	else
		samples=$(sed '/^#/d' ${inputs["meta"]} | awk '{print $3}')
		for i in ${samples[@]}; do
			echo "-V ${vcf_dir}/${inputs['cohort_id']}.$i.${input_gvcf_stage}.g.vcf"
		done > ${inputs['tmp_prefix']}gvcf.list
	fi

	[[ -s ${inputs['tmp_prefix']}gvcf.list ]] || { echo "\n\e[38;5;1mERROR\e[0m: no gvcf files were found, please check and try again"; return 1; }

}