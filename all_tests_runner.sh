#!/bin/bash
rm simrun
rm hc_*
rm mass_???
rm ms??
rm *.ps
rm *.pdf

Rscript hard_case_CV_airbnb.R &> hc_Ca &
Rscript hard_case_DMRnet_insurance.R &> hc_Di &
Rscript hard_case_DMRnet_simulations.R &> hc_Ds &
Rscript hard_case_GLAMER_promoter.R &> hc_Gp &
Rscript hard_case_LAPACK_SVD_insurance.R &> hc_LSi &
Rscript hard_case_NA_insurance.R &> hc_Ni &
Rscript hard_case_permute_columns.R &> hc_pc &
Rscript hard_case_SOSnet.R &> hc_S &

Rscript massive_adult.R &> mass_adu &
Rscript massive_antigua.R &> mass_ant &
Rscript massive_insurance.R &> mass_ins &
Rscript massive_promoter.R &> mass_pro &

./simmulation_runner.sh &> simrun &
