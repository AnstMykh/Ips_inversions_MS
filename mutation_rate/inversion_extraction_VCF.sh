!/bin/sh
cd $1   

list_bed=$(ls *.bed)

#MJ, MN and R are the names of haplotypes depending on their frequency: MJ - major, MN -minor, R - 'recombinants' (they were not included in the frequency count)  
for i in $list_bed
    do
    name=$(echo $i| cut -d '.' -f1)
    echo $name
    MJ=$(ls *_MJ.args | grep $name)
    echo $MJ
    MN=$(ls *_MN.args | grep $name)
    echo $MN
    R=$(ls *_R.args | grep $name)
    echo $R

    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/Contig_VCFs/new/Bark_beetle_merged_filtered_no_outliers_fin.vcf.gz \
    -L $i \
    -sn $1/$MJ \
    -O ~/mut_load/OUTPUT/outliers/$name'_MJ_tmp.vcf'

 #--exclude-non-variants \
    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/mut_load/OUTPUT/outliers/$name'_MJ_tmp.vcf' \
    --remove-unused-alternates \
    -O ~/mut_load/OUTPUT/outliers/$name'_MJ.vcf'



    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/Contig_VCFs/new/Bark_beetle_merged_filtered_no_outliers_fin.vcf.gz \
    -L $i \
    -sn $1/$MN \
    -O ~/mut_load/OUTPUT/outliers/$name'_MN_tmp.vcf'

 #--exclude-non-variants \
    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/mut_load/OUTPUT/outliers/$name'_MN_tmp.vcf' \
    --remove-unused-alternates \
    -O ~/mut_load/OUTPUT/outliers/$name'_MN.vcf'



    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/Contig_VCFs/new/Bark_beetle_merged_filtered_no_outliers_fin.vcf.gz \
    -L $i \
    -sn $1/$R \
    -O ~/mut_load/OUTPUT/outliers/$name'_R_tmp.vcf'


 #--exclude-non-variants \
    gatk SelectVariants \
    -R ~/Ref_fas/SortedAndRenamed.clean.fasta \
    -V ~/mut_load/OUTPUT/outliers/$name'_R_tmp.vcf' \
    --remove-unused-alternates \
    -O ~/mut_load/OUTPUT/outliers/$name'_R.vcf'
    done 
