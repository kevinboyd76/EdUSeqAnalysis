# 1 - Generate Sigma Values for Adjust and Sample files
sbatch eduSeq_sigmaCalculation.sh ../results/aligned/EduHU_HCT_TotalSheared_set2A.sam ../results/aligned/EduHU_HCT_Biotin_set2A.sam

# 2 - Generate Sigma Output and Sigma Quality Report
sbatch sigmaCalc_submit.sh

# 3 - Normalize Sigma
sbatch sigmaNorm.sh
