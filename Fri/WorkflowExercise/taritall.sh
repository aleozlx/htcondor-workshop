#!/bin/sh

/bin/tar -cvzf results_$1.tar.gz perm_combined_$1.Rdat perm_summary_$1.txt sim_geno_results_$1.Rdat find_model_results_$1.Rdat qtl_$1.Rdat refined_qtl_$1.Rdat refined_qtl_summary_$1.txt fit_qtl_results_$1.Rdat  > tarall_$1.log
