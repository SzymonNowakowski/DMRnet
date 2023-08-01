#!/bin/bash



echo "============================================================" > tag_text
echo "hard_case_CV_airbnb.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_Ca >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_DMRnet_insurance.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_Di >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_DMRnet_simulations.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_Ds >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_GLAMER_promoter.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_Gp >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_LAPACK_SVD_insurance.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_LSi >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_NA_insurance.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_Ni >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_permute_columns.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_pc >> tag_text

echo "============================================================" >> tag_text
echo "hard_case_SOSnet.R" >> tag_text
echo "============================================================" >> tag_text
cat hc_S >> tag_text


