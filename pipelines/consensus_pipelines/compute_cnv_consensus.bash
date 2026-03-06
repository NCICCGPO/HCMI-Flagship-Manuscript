#1/bin/bash

#Set things up
#requires R-4.1 in path

#Store important information as variables
cd ${OUTPUT_DIR}
sample_name=$(tail -n +2 absolute_solutions_for_forcecalling.txt | head -n $SGE_TASK_ID | tail -n 1 | awk '{print $2}')
absolute_solution=$(tail -n +2 absolute_solutions_for_forcecalling.txt | head -n $SGE_TASK_ID | tail -n 1 | awk '{print $8}')



#Forcecall absolute
cd ${ABSOLUTE_OUTPUT_PATH}/$sample_name
Rscript ABSOLUTE_extract_cli_start.R --solution_num $absolute_solution --analyst_id force_called --rdata_modes_fn $sample_name.PP-modes.data.RData --sample_name $sample_name --results_dir . --abs_lib_dir ${ABSOLUTE_SOFTWARE_PATH}



#Move the seg to the output directory
cp ${ABSOLUTE_OUTPUT_PATH}/$sample_name/reviewed/SEG_MAF/$sample_name.segtab.txt ${ABSOLUTE_OUTPUT_PATH}