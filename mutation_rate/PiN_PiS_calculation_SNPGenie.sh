####reverse strand data preparation 


### FASTA
for i in {1..36}
  do 
  perl fasta2revcom.pl IpsContig$i.fa 
  done >> log.txt


#### GTF

while IFS= read -r line
    do
    contig=$(echo $line | cut -d " " -f 1)
    length=$(echo $line | cut -d " " -f 2)
    perl gtf2revcom.pl gtf_$contig.gtf $length
    done < fasta_lengths.txt




###VCF
## outlying regions
cd ~/mut_load/OUTPUT/outliers/tmp                           

list=$(ls *.vcf)

for i in $list
    do
    contig=$(echo $i | cut -d "_" -f 1)
    length=$(grep $contig ~/piNpiS/fasta_lengths.txt  | cut -d " " -f 2)
    perl ~/piNpiS/vcf2revcom.pl ~/mut_load/OUTPUT/outliers/tmp/$i $length
    done


## collinear

cd ~/mut_load/OUTPUT/collinear_sample/ 
list=$(ls *_collinear.vcf)

for i in $list
    do
    contig=$(echo $i | cut -d "_" -f 1)
    length=$(grep $contig" " ~/piNpiS/fasta_lengths.txt | cut  -d " " -f 2)
    perl ~/piNpiS/vcf2revcom.pl ~/mut_load/OUTPUT/collinear_sample/$i $length
    done 




#################### MUTATION LOAD #############################


####whole-genome mut_load calculations: collinear 

# "+" strand
cd ~/mut_load/OUTPUT/collinear_sample/
list=$(ls *_collinear.vcf)
cd ~/piNpiS
for i in $list
  do 
  contig=$(echo $i | cut -d "_" -f 1)
  perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/collinear_sample/$i --fastafile=$contig.fa --gtffile=gtf_$contig.gtf --outdir $contig'_collinear'
  #perl ~/piNpiS/snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/collinear_sample/$contig'_collinear_revcom'.vcf --fastafile=~/piNpiS/$contig'_revcom'.fa --gtffile=~/piNpiS/gtf_$contig'_revcom'.gtf --outdir $contig'_revcom_collinear'
  done


# "-" strand
cd ~/mut_load/OUTPUT/collinear_sample/
list=$(ls *_revcom_collinear.vcf)
cd ~/piNpiS
for i in $list
  do 
  contig=$(echo $i | cut -d "_" -f 1)
  #perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/collinear_sample/$i --fastafile=$contig.fa --gtffile=gtf_$contig.gtf --outdir $contig'_collinear'
  perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/collinear_sample/$i --fastafile=$contig'_revcom'.fasta --gtffile=gtf_$contig'_revcom'.gtf --outdir $contig'_revcom_collinear'
  done




####whole-genome mut_load calculations: inversions 
cd ~/mut_load/OUTPUT/outliers/

list=$(ls *.vcf)

# "+" strand
cd ~/piNpiS
for i in $list
  do
  contig=$(echo $i | cut -d "_" -f 1)
  name=$(echo $i | cut -d "." -f 1)
  perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/outliers/$i --fastafile=$contig.fa --gtffile=gtf_$contig.gtf --outdir $name
  #perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/collinear_sample/$contig'_collinear_revcom'.vcf --fastafile=$contig'_revcom'.fasta --gtffile=gtf_$contig'_revcom'.gtf --outdir $contig'_revcom_collinear'

  done

#### "-" strand 
cd ~/mut_load/OUTPUT/outliers/revcom/

list=$(ls *.vcf)
cd ~/piNpiS

for i in $list
  do 
  contig=$(echo $i | cut -d "_" -f 1)
  name=$(echo $i | cut -d "." -f 1)
  #perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/outliers/$i --fastafile=$contig.fa --gtffile=gtf_$contig.gtf --outdir $name
  perl snpgenie.pl --vcfformat=1 --snpreport=/home/anastasiia.mykhailenko/mut_load/OUTPUT/outliers/revcom/$i --fastafile=$contig'_revcom'.fasta --gtffile=gtf_$contig'_revcom'.gtf --outdir $name
  done

