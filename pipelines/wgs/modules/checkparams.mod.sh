#
#  checkparams (module)
#
#	summary:
#   + check_input_json()
#   + check_id()
#	+ check_ref()
#	+ check_bwa_idx()
#	+ check_gatk_dict()
#	+ check_samtools_fai()
#
check_id() {
    local \
        id_type=$1 \
        id=$2 \
        status=0

    # check the id is not empty or null
    [[ -z $id ]] && { echo 'ERROR: input $id_type id is empty, please choose another'; status=1; }

    # check the id type has the relevant identifier: c_..., s_...
    if [[ $id_type == 'cohort' ]]; then
        [[ $id == "c_"* ]] || { echo "ERROR: cohort ids must start with c_... please choose another"; status=1; }
    elif [[ $id_type == "s_"* ]]; then
        [[ $id == "s_"* ]] || { echo "ERROR: sample ids must start with s_... please choose another"; status=1; }
    fi

    # check the id has no spaces or tabs

    return $status
}

function check_ref() {
       if [[ "${inputs["ref"]}" == NULL ]]; then
          echo -e " -r,--ref not provided! Exiting..."; 1>&2;
          return 1
       elif [ ! -f "${inputs["ref"]}" -o ! -s "${inputs["ref"]}" ]; then
          echo -e " Problem with reference file. Please check that it exists and is not empty..."; 1>&2;
          return 1
       fi
}
function check_bwa_idx() {
  if [[ ! -f "${ref_dir}/${inputs['ref_label']}/bwa.${inputs['ref_label']}.amb" ]]; then
    echo " can't find the bwa index files for ${inputs["ref_label"]}"
		return 1
  fi
}

function check_gatk_dict() {
	if [[ ! -f "${ref_dir}/${inputs['ref_label']}/${inputs['ref_label']}.dict" ]]; then
		echo " can't find gatk reference dictionary"
		return 1
	fi
}
function check_samtools_fai() {
       if [[ ! -f "${ref_dir}/${inputs['ref_label']}/${inputs['ref_label']}.fai" ]]; then
          echo " can't find samtools fai index"
          return 1
       fi
}