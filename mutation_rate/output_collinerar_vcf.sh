!/bin/sh
cd ~/Inversions/PCA_outliers/wide_coordinates/

list_bed=$(ls *.bed)

for i in $list_bed
    do
    contig=$(echo $i | cut -d "_" -f1 )

    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/Contig_VCFs/new/GATK_VCF_Ips_contig/$contig.vcf \
    -XL $i \
    -O ~/mut_load/OUTPUT/$contig'_collinear_tmp'.vcf

    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/mut_load/OUTPUT/$contig'_collinear_tmp'.vcf \
    --remove-unused-alternates \
    -O ~/mut_load/OUTPUT/$contig'_collinear'.vcf
    done 
