#$_1: r150_bed_A $_2: r150_bed_B $_3: bedpe_A $_4: bedpe_B $_5: Loop_A $_6: Loop_B $_7: FDR_thr

workdir=`pwd`

cd `dirname $0`
scriptdir=`pwd`

cd $workdir

echo "${scriptdir}/merge_tracloop merged_loops.txt merged_anchors.txt $6 $7"
${scriptdir}/merge_tracloop merged_loops.txt merged_anchors.txt $6 $7

echo "${scriptdir}/count_tag_on_region merged_anchors.txt $1 Acc_on_merged_anchors_A.txt"
${scriptdir}/count_tag_on_region	merged_anchors.txt $1 Acc_on_merged_anchors_A.txt

echo "${scriptdir}/count_tag_on_region merged_anchors.txt $2 Acc_on_merged_anchors_B.txt"
${scriptdir}/count_tag_on_region merged_anchors.txt $2 Acc_on_merged_anchors_B.txt

echo "${scriptdir}/get_petc_in_loop_from_bedpe merged_loops.txt $3 PETc_on_merged_loops_A.txt"
${scriptdir}/get_petc_in_loop_from_bedpe merged_loops.txt $3 PETc_on_merged_loops_A.txt

echo "${scriptdir}/get_petc_in_loop_from_bedpe merged_loops.txt $4 PETc_on_merged_loops_B.txt"
${scriptdir}/get_petc_in_loop_from_bedpe merged_loops.txt $4 PETc_on_merged_loops_B.txt

echo "${scriptdir}/call_diff_loop PETc_on_merged_loops_A.txt PETc_on_merged_loops_B.txt Acc_on_merged_anchors_A.txt Acc_on_merged_anchors_B.txt merged_loops_comp_summary"
${scriptdir}/call_diff_loop PETc_on_merged_loops_A.txt PETc_on_merged_loops_B.txt Acc_on_merged_anchors_A.txt Acc_on_merged_anchors_B.txt merged_loops_comp_summary

echo "python ${scriptdir}/cal_qvalue_diff_tracloop.py merged_loops_comp_summary merged_loops_comp_summary2"
python ${scriptdir}/cal_qvalue_diff_tracloop.py merged_loops_comp_summary merged_loops_comp_summary2

echo "python ${scriptdir}/get_significant_diff_tracloop.py merged_loops_comp_summary2 ${FDR_thr} merged_loops_increased_FDR${FDR_thr} merged_loops_decreased_FDR${FDR_thr} merged_loops_nochange_FDR${FDR_thr}"
python ${scriptdir}/get_significant_diff_tracloop.py merged_loops_comp_summary2 ${FDR_thr} merged_loops_increased_FDR${FDR_thr} merged_loops_decreased_FDR${FDR_thr} merged_loops_nochange_FDR${FDR_thr} 
