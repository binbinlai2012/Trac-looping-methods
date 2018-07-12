#$_1: bedpe  $_2: peak_file(MACS/SICER)  $_3: tag  $_4: chrlen  $_5: win  $_6: minLength  $_7: maxLength  $_8: min_PETs  $_9: FDR_thr

# call interaction background model
${PWD}/call_interaction_model -L $4 -b $1 -p $3 -w $5 -l $6 -u $7 -R $2 -d 0.1

Rscript ${PWD}/smooth_spline.r ${3}.distance.dist.txt 0.25 ${3}.distance.dist.sm.txt

Rscript ${PWD}/smooth_spline.r ${3}.distance.PETs_dist.txt 0.5 ${3}.distance.PETs_dist.sm.txt

Rscript ${PWD}/smooth_spline.r ${3}.accm.dist.txt 0.25 ${3}.accm.dist.sm.txt 

Rscript ${PWD}/smooth_spline.r ${3}.accm.PETs_dist.txt 0.25 ${3}.accm.PETs_dist.sm.txt 

# call interaction
${PWD}/call_interaction -R $2 -b $1 -p $3 -c 1 -w $5 -l $6 -u $7 \
	-0 ${3}.distance.dist0.txt \
	-1 ${3}.distance.dist.sm.txt \
	-2 ${3}.distance.PETs_dist0.txt \
	-3 ${3}.distance.PETs_dist.sm.txt \
	-4 ${3}.accm.dist.sm.txt \
	-5 ${3}.accm.PETs_dist.sm.txt \
	-6 ${3}.total_valid_space_size.txt \
	-7 ${3}.total_valid_PETs.txt 
#	-8 ${3}.total_valid_PETs_in_anchor.txt 

# call FDR
python ${PWD}/cal_qvalue.py ${3}.interaction.summary ${3}.interaction.summary.fdr
python ${PWD}/flt_sig_interaction.py ${3}.interaction.summary.fdr ${8} ${9} ${3}.interaction.c${8}.fdr${9}

# merge significant interaction
${PWD}/merge_sig_interaction ${3}.interaction.c${8}.fdr${9} $1 $6 $7 ${3}.total_valid_PETs.txt $3
