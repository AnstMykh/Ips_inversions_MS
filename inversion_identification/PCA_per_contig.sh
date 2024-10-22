#parcing Ips genome into separape contigs .vcf

tabix --list-chroms whole_genome.vcf.gz > IpsContigs.txt
while IFS= read -r line; do   
    vcftools --gzvcf whole_genome.vcf.gz  --chr $line --recode --recode-INFO-all --out $line; # you might want to rename your files later
    done < IpsContigs.txt



# performing PCA on each contig, move all .vcf to the empty folder:
# making a list of all files
list=$(ls .vcf)

#looping through the list of contigs
for i in $list:
  do
  ind=$(echo $i |cut -d"." -f 1)
  plink --vcf $i --pca header --allow-extra-chr --out $ind  # optional: --mind 0.5
  done
