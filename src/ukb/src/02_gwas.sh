## Null models

run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ; \
    ls *.bed | sed -e 's/.bed//g' > files_to_merge.txt; \
    plink --merge-list files_to_merge.txt --make-bed --autosome-xy --out \
    ukb22418_c1_22_v2_merged; \
    rm files_to_merge.txt;"

dx run swiss-army-knife \
    -iin="/GWAS/main/pheno_brain_pcs.tsv" \
    -icmd="${run_merge}" --tag="Step1" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="/GWAS/results_brain_pcs" --brief --yes


#############
data_file_dir="/GWAS/results_brain_pcs"
pheno_dir="/GWAS/main"

run_plink_qc="plink2 --bfile ukb22418_c1_22_v2_merged \
--keep pheno_brain_pcs.tsv --autosome \
--maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 \
--mind 0.1 --write-snplist --write-samples \
--no-id-header --out WES_array_snps_qc_pass"

dx run swiss-army-knife \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam" \
    -iin="${pheno_dir}/pheno_brain_pcs.tsv" \
    -icmd="${run_plink_qc}" \
    --tag="Step1" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${data_file_dir}/" \
    --brief --yes

#############
data_file_dir="/GWAS/results_brain_pcs"
pheno_dir="/GWAS/main"

## --phenoCol Biotype \
## --bt

run_regenie_step1="regenie --step 1 \
--lowmem --out tinnitus_results --bed ukb22418_c1_22_v2_merged \
--phenoFile pheno_brain_pcs.tsv \
--covarFile covariates_brainPC_with_srt.tsv \
--extract WES_array_snps_qc_pass.snplist \
--covarCol age --covarCol genotype --covarCol srt --covarCol eTIV \
--covarCol PC_1 --covarCol PC_2 --covarCol PC_3 --covarCol PC_4 --covarCol PC_5 \
--covarCol PC_6 --covarCol PC_7 --covarCol PC_8 --covarCol PC_9 --covarCol PC_10 \
--bsize 1000 --loocv --gz --threads 16"

# --- EXECUTION ---
dx run swiss-army-knife \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bed" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.bim" \
    -iin="${data_file_dir}/ukb22418_c1_22_v2_merged.fam" \
    -iin="${data_file_dir}/WES_array_snps_qc_pass.snplist" \
    -iin="${pheno_dir}/pheno_brain_pcs.tsv" \
    -iin="${pheno_dir}/covariates_brainPC_with_srt.tsv" \
    -icmd="${run_regenie_step1}" \
    --tag="Regenie_Step3_Tinnitus" \
    --name "Tinnitus_Step3_Model" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${data_file_dir}/" \
    --brief --yes

#############
exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - interim 450k release/"
data_file_dir="/GWAS/results_brain_pcs"
pheno_dir="/GWAS/main"
field_name="ukb23149"

for chr in {1..22}; do
    run_plink_wes="plink2 --bfile ${field_name}_c${chr}_b0_v1 \
        --no-pheno --keep pheno_brain_pcs.tsv --autosome \
        --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-15 --mind 0.1 \
        --write-snplist --write-samples --no-id-header \
        --out WES_c${chr}_snps_qc_pass"

        dx run swiss-army-knife \
        -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bed" \
        -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bim" \
        -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.fam" \
        -iin="${pheno_dir}/pheno_brain_pcs.tsv" \
        -icmd="${run_plink_wes}" \
        --tag="Step2_QC" \
        --name "QC_WES_chr${chr}" \
        --instance-type "mem1_ssd1_v2_x16" \
        --destination="${data_file_dir}/" \
        --brief --yes
done


#############
exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - interim 450k release/"
data_file_dir="/GWAS/results_brain_pcs"
pheno_dir="/GWAS/main"
field_name="ukb23149"

for chr in {1..22}; do
    run_regenie_cmd="regenie --step 2 --bed ${field_name}_c${chr}_b0_v1 --out tinnitus_assoc.c${chr} \
        --phenoFile pheno_subtype.tsv --covarFile covariates_subtype_with_srt.tsv \
        --bt --approx --firth-se --firth --extract WES_c${chr}_snps_qc_pass.snplist \
        --phenoCol Biotype \
        --covarCol age --covarCol genotype --covarCol srt \
        --covarCol PC_1 --covarCol PC_2 --covarCol PC_3 --covarCol PC_4 --covarCol PC_5 \
        --covarCol PC_6 --covarCol PC_7 --covarCol PC_8 --covarCol PC_9 --covarCol PC_10 \
        --pred tinnitus_results_pred.list --bsize 200 \
        --pThresh 0.05 --minMAC 3 --threads 16 --gz"

        dx run swiss-army-knife \
        -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bed" \
        -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bim" \
        -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.fam" \
        -iin="${data_file_dir}/WES_c${chr}_snps_qc_pass.snplist" \
        -iin="${pheno_dir}/pheno_subtype.tsv" \
        -iin="${pheno_dir}/covariates_subtype_with_srt.tsv" \
        -iin="${data_file_dir}/tinnitus_results_pred.list" \
        -iin="${data_file_dir}/tinnitus_results_1.loco.gz" \
        -icmd="${run_regenie_cmd}" --tag="Step2_Assoc" --instance-type "mem1_ssd1_v2_x16" \
        --name "regenie_step2_chr${chr}" \
        --destination="${data_file_dir}/" --brief --yes
done

#############
merge_cmd='out_file="tinnitus_assoc_merged.txt"
  cp /mnt/project/GWAS/results_subtype/*.regenie.gz .
  gunzip *.regenie.gz
  echo -e "CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA" > $out_file
  files="./*.regenie"
for f in $files
do
    tail -n+2 $f | tr " " "\t" >> $out_file
done

rm *.regenie'

# Define your data directory
data_file_dir="/GWAS/results_subtype"

# Execute the merge
dx run swiss-army-knife \
    -iin="${data_file_dir}/tinnitus_assoc.c1_Biotype.regenie.gz" \
    -icmd="${merge_cmd}" \
    --tag="Final_Merge" \
    --instance-type "mem1_ssd1_v2_x16" \
    --destination="${data_file_dir}" \
    --brief --yes

############# only for brain volumes PCs
exome_file_dir="/Bulk/Exome sequences/Population level exome OQFE variants, PLINK format - interim 450k release/"
data_file_dir="/GWAS/results_brain_pcs"
pheno_dir="/GWAS/main"
field_name="ukb23149"

for pheno in PC1_brain PC2_brain; do
    for chr in {1..22}; do
        run_regenie_cmd="regenie --step 2 \
            --bed ${field_name}_c${chr}_b0_v1 \
            --out ${pheno}_assoc.c${chr} \
            --phenoFile pheno_brain_pcs.tsv \
            --covarFile covariates_brainPC_with_srt.tsv \
            --extract WES_c${chr}_snps_qc_pass.snplist \
            --phenoCol ${pheno} \
            --covarCol age --covarCol genotype --covarCol srt --covarCol eTIV \
            --covarCol PC_1 --covarCol PC_2 --covarCol PC_3 --covarCol PC_4 --covarCol PC_5 \
            --covarCol PC_6 --covarCol PC_7 --covarCol PC_8 --covarCol PC_9 --covarCol PC_10 \
            --pred tinnitus_results_pred.list \
            --bsize 200 \
            --minMAC 3 \
            --threads 16 \
            --gz"

        dx run swiss-army-knife \
            -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bed" \
            -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.bim" \
            -iin="${exome_file_dir}/${field_name}_c${chr}_b0_v1.fam" \
            -iin="${data_file_dir}/WES_c${chr}_snps_qc_pass.snplist" \
            -iin="${pheno_dir}/pheno_brain_pcs.tsv" \
            -iin="${pheno_dir}/covariates_brainPC_with_srt.tsv" \
            -iin="${data_file_dir}/tinnitus_results_pred.list" \
            -iin="${data_file_dir}/tinnitus_results_1.loco.gz" \
            -iin="${data_file_dir}/tinnitus_results_2.loco.gz" \
            -icmd="${run_regenie_cmd}" \
            --tag="Step2_Assoc" \
            --instance-type "mem1_ssd1_v2_x16" \
            --name "regenie_${pheno}_chr${chr}" \
            --destination="${data_file_dir}/" \
            --brief --yes
    done
done

############# only for brain PC
for pheno in PC1_brain PC2_brain; do

    merge_cmd="out_file=\"${pheno}_assoc_merged.txt\"
    cp /mnt/project/GWAS/results_brain_pcs/${pheno}_assoc.c*.regenie.gz .
    gunzip *.regenie.gz
    echo -e \"CHROM\tGENPOS\tID\tALLELE0\tALLELE1\tA1FREQ\tN\tTEST\tBETA\tSE\tCHISQ\tLOG10P\tEXTRA\" > \$out_file
    files=\"./*.regenie\"
    for f in \$files
    do
        tail -n+2 \$f | tr \" \" \"\t\" >> \$out_file
    done

    rm *.regenie"

    # Define your data directory
    data_file_dir="/GWAS/results_brain_pcs"

    # Execute the merge
    dx run swiss-army-knife \
        -iin="${data_file_dir}/${pheno}_assoc.c1_${pheno}.regenie.gz" \
        -icmd="${merge_cmd}" \
        --tag="Final_Merge" \
        --name "merge_${pheno}" \
        --instance-type "mem1_ssd1_v2_x16" \
        --destination="${data_file_dir}" \
        --brief --yes

done

