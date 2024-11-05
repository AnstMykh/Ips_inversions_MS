list=$(ls *.vcf)

for i in $list
    do
    name=$(echo $i | cut -d "." -f 1)
    vcftools --vcf path/to/$i --thin 10000 --maf 0.05 --recode --out $name
    grep -v "^##" $name.recode.vcf | cut -f2 > $name.pos.txt
    grep -v "^##" $i | cut -f2 > path/to/$name.pos.txt
    plink --vcf $i  --make-bed  --allow-extra-chr --out ./new/$name
    plink --bfile $name --cow --r2 square --allow-extra-chr --out path/to/$name
    done 
