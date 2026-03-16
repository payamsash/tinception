## Null models

run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ; \
    ls *.bed | sed -e 's/.bed//g' > files_to_merge.txt; \
    plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out \
    ukb22418_c1_22_v2_merged; \
    rm files_to_merge.txt;"

dx run swiss-army-knife \
    -iin="/GWAS/main/pheno_case_control.tsv" \
    -icmd="${run_merge}" --tag="Step1" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="/GWAS/results" --brief --yes


#############
# Define the project-relative paths
data_file_dir="/GWAS/results"
pheno_dir="/GWAS/main"

run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged \
--keep pheno_case_control.tsv --autosome \
--maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 \
--mind 0.1 --write-snplist --write-samples \
--no-id-header --out WES_array_snps_qc_pass"

dx run swiss-army-knife \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam" \
    -iin="${pheno_dir}/pheno_case_control.tsv" \
    -icmd="${run_plink_qc}" \
    --tag="Step1" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${data_file_dir}/" \
    --brief --yes

#############
data_file_dir="/GWAS/results"
pheno_dir="/GWAS/main"

run_regenie_step1="regenie --step 1 \
--lowmem --out tinnitus_results --bed ukb22418_c1_22_v2_merged \
--phenoFile pheno_case_control.tsv \
--covarFile covariates_case_control_with_srt.tsv \
--extract WES_array_snps_qc_pass.snplist \
--phenoCol tin_status \
--covarCol age --covarCol genotype --covarCol srt \
--covarCol PC_1 --covarCol PC_2 --covarCol PC_3 --covarCol PC_4 --covarCol PC_5 \
--covarCol PC_6 --covarCol PC_7 --covarCol PC_8 --covarCol PC_9 --covarCol PC_10 \
--bsize 1000 --bt --loocv --gz --threads 16"

# --- EXECUTION ---
dx run swiss-army-knife \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam" \
    -iin="${data_file_dir}/WES_array_snps_qc_pass.snplist" \
    -iin="${pheno_dir}/pheno_case_control.tsv" \
    -iin="${pheno_dir}/covariates_case_control_with_srt.tsv" \
    -icmd="${run_regenie_step1}" \
    --tag="Regenie_Step1_Tinnitus" \
    --name "Tinnitus_Step1_Model" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${data_file_dir}/" \
    --brief --yes

#############
