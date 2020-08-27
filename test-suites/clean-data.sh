#
#  CLEAN DATA - maintenence script for pipeline use
#
#	* clean all intermediate and final TEST SUITE output data
#
#	Jack Morrice
#
###############################################################

#!/bin/bash
cd '../'

source includes/locations.sh
source includes/utilities.sh

workflow() { 
	local argv=("$@")

	./sap.sh utils clean-data

	[[ -d ${rds_dir} ]] && echo '  cleaning reads files' && rm ${rds_dir}/* 2> /dev/null
	[[ -d ${rds_dir} ]] && echo '  cleaning truth vcf files' && rm ${vcf_dir}/*.truth.* 2> /dev/null
}

#
#  tasks
#

#
#  run workflow
#
workflow "$@"