#1.Generation of check sums. Run in parent directory to the directories where sequencing data are stored. Compare results with data provided by the sequencing facility.#
for folder in */
	do
	echo $folder
	cd $folder
	for file in *.fastq.gz
	do
	echo $file
	md5sum $file >> ~/TK-2760/Md5sums_generated.txt
	done
	cd ..
done

#2.Quality control with FastQC. Run in directory where the sequencing data are#
#2.1.Quality control in the whole data at once ie. all files taken together - to assess at once what is the overall data quality#

for folder in */
	do
	echo $folder
	cd $folder
	zcat *.fastq.gz >> ~/TK-2760/All_data_temporal.fq
	cd ..
done

cd ..

cat All_data_temporal.fq | ~/FastQC/fastqc stdin -t 40

rm All_data_temporal.fq

#2.2.Per sample quality control - to check if there is variation between samples, to see which sampels performed poorly#

for folder in */
	do
	echo $folder
	cd $folder
	zcat *.fastq.gz | ~/FastQC/fastqc stdin -t 10
	cd ..
done

#3.Read number calculation - to learn how many reads we have per individual and to check if numer of R1 & R2 reads is the same within individual - should be!#

for folder in */
	do
	echo $folder
	cd $folder
	for file in *_R1_001.fastq.gz
	do
	echo $file >> ~/TK-2760/R1_reads.txt
	zcat $file | echo $((`wc -l`/4)) >> ~/TK-2760/R1_reads.txt
	done
	for file in *_R2_001.fastq.gz
	do
	echo $file >> ~/TK-2760/R2_reads.txt
	zcat $file | echo $((`wc -l`/4)) >> ~/TK-2760/R2_reads.txt
	done
	cd ..
done

#4.Copy files to one common directory AND chnge their names#

for folder in */
	do
	echo $folder
	cd $folder
	for file in *_R1_001.fastq.gz
	do
	nameis=$(echo $folder | cut -d'-' -f 3-6 | tr -d '/')
	echo $nameis
	zcat $file >> ~/korniki/$nameis"_R1_001.fastq"
	done
	cd ..
done

for folder in */
	do
	echo $folder
	cd $folder
	for file in *_R2_001.fastq.gz
	do
	nameis=$(echo $folder | cut -d'-' -f 3-6 | tr -d '/')
	echo $nameis
	zcat $file >> ~/korniki/$nameis"_R2_001.fastq"
	done
	cd ..
done

#5.Trimming adaptors & low quality reads with Trimmomatic#

lista=$(ls *.fastq | cut -d"_" -f 1-3 --output-delimiter="_" | sort | uniq)

for indiv in $lista
	do
	echo "mapping "$indiv
	ind=$(echo $indiv |cut -d"_" -f 1)
	echo $ind
	java -Xmx32g -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 $ind"_R1_001.fastq" $ind"_R2_001.fastq" $ind"_R1_001.fq.gz" $ind"_forward_unpaired.fq.gz" $ind"_R2_001.fq.gz"  $ind"_reverse_unpaired.fq.gz" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
	done

#6.Read number calculation - to learn how many reads we lost while trimming#

for file in *R1_001.fq.gz
	do
	echo $file >> R1_reads_after_trim.txt
	zcat $file | echo $((`wc -l`/4)) >> R1_reads_after_trim.txt
	done
for file in *R2_001.fq.gz
	do
	echo $file >> R2_reads_after_trim.txt
	zcat $file | echo $((`wc -l`/4)) >> R2_reads_after_trim.txt
	done

#7.Maping with bowtie2 AND .sam to .bam conversion with samtools#
bowtie2-build SortedAndRenamed.clean.fasta ref_idx

lista=$(ls *R1_001.fq.gz | cut -d"_" -f 1-3 --output-delimiter="_" | sort | uniq)

for indiv in $lista
	do
	echo "mapping "$indiv
	ind=$(echo $indiv |cut -d"_" -f 1)
	echo $ind
	bowtie2 -x ref_idx -1 $ind"_R1_001.fq.gz" -2 $ind"_R2_001.fq.gz" -S $ind".sam" --rg-id $ind --rg SM:$ind -p 90 --no-unal 2>> Bowtie_Mapping_Report.txt
	samtools view -bS $ind".sam" > $ind".bam"
	rm $ind".sam"
	done

#8.Refernce indexing & dictionary creation for downstream analysis with samtools and picard#
samtools faidx SortedAndRenamed.clean.fasta
picard CreateSequenceDictionary R=SortedAndRenamed.clean.fasta O=SortedAndRenamed.clean.dict

#9. Bam file sorting & indexing with samtools#
files=*.bam

for plik in $files
	do
	nameis=$(echo $plik | cut -d'.' -f 1)
	echo $nameis
	samtools sort -T tmp -o $nameis"_coordsorted.bam" $nameis".bam"
	samtools index $nameis"_coordsorted.bam"
	done

#10. Mapping quality check with Qualimap#
mkdir QualimapResults

for file in *_coordsorted.bam
	do
	nameis=$(echo $file | cut -d'_' -f 1)
	unset DISPLAY
	qualimap bamqc -bam $file -nt 80 --java-mem-size=80G -outdir QualimapResults -outfile $nameis.pdf -outformat PDF
	done

#11. Duplicate read removal#

for file in *_coordsorted.bam
	do
	nameis=$(echo $file | cut -d'_' -f 1)
	picard MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=$file OUTPUT=$nameis"_coordsorted_rmdups.bam" METRICS_FILE=$nameis"_picard_markduplicates_metrics.txt" REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true
	done

#12. Calculate reads in bams that have removed duplicates#
for file in *_coordsorted_rmdups.bam
	do
	echo $file
	samtools idxstats $file | awk -F '\t' '{sum+=$3}END{print sum}'
	done > Mapped_reads_rmdups.txt

#13. Variant calling was first performed to create database of known SNPs to use it further with BQSR and VQSR#
#13.1. Emiting single sample GVCFs with GATK Haplotype caller#

for file in *_coordsorted_rmdups.bam
	do
	echo $file
	nameis=$(echo $file | cut -d'_' -f 1)
	echo $nameis
	gatk HaplotypeCaller \
	-R SortedAndRenamed.clean.fasta \
	-I $file \
	-O $nameis.g.vcf.gz \
	--heterozygosity 0.01 \
	-ERC GVCF
	done

#13.2. Combine variants was run after HaplotypeCaller and before GenotypeGVCFs#
gatk CombineGVCFs \
 -R SortedAndRenamed.clean.fasta \
 --variant AAS10.g.vcf.gz \
 --variant AAS18.g.vcf.gz \
 --variant AAS2.g.vcf.gz \
 --variant AAS21.g.vcf.gz \
 --variant AAS22.g.vcf.gz \
 --variant AAS24.g.vcf.gz \
 --variant AAS25-DUP.g.vcf.gz \
 --variant AAS25.g.vcf.gz \
 --variant AAS28.g.vcf.gz \
 --variant AAS29.g.vcf.gz \
 --variant AAS3.g.vcf.gz \
 --variant AAS30.g.vcf.gz \
 --variant AAS41.g.vcf.gz \
 --variant AAS6.g.vcf.gz \
 --variant ASA10.g.vcf.gz \
 --variant ASA11.g.vcf.gz \
 --variant ASA12.g.vcf.gz \
 --variant ASA124.g.vcf.gz \
 --variant ASA125.g.vcf.gz \
 --variant ASA126.g.vcf.gz \
 --variant ASA13.g.vcf.gz \
 --variant ASA14.g.vcf.gz \
 --variant ASA2.g.vcf.gz \
 --variant ASA20.g.vcf.gz \
 --variant ASA21.g.vcf.gz \
 --variant ASA46.g.vcf.gz \
 --variant ASA47.g.vcf.gz \
 --variant ASA8.g.vcf.gz \
 --variant BAD1.g.vcf.gz \
 --variant BAD16.g.vcf.gz \
 --variant BAD17.g.vcf.gz \
 --variant BAD18.g.vcf.gz \
 --variant BAD3.g.vcf.gz \
 --variant BAW148.g.vcf.gz \
 --variant BAW150.g.vcf.gz \
 --variant BAW152.g.vcf.gz \
 --variant BAW153.g.vcf.gz \
 --variant BAW154.g.vcf.gz \
 --variant BAW156.g.vcf.gz \
 --variant BAW158.g.vcf.gz \
 --variant BAW160.g.vcf.gz \
 --variant BAW162.g.vcf.gz \
 --variant BAW164.g.vcf.gz \
 --variant BAW165.g.vcf.gz \
 --variant BAW175.g.vcf.gz \
 --variant BAW178.g.vcf.gz \
 --variant BIL1-DUP.g.vcf.gz \
 --variant BIL10.g.vcf.gz \
 --variant BIL11.g.vcf.gz \
 --variant BIL20.g.vcf.gz \
 --variant BIL27.g.vcf.gz \
 --variant BIL34.g.vcf.gz \
 --variant BIL40.g.vcf.gz \
 --variant BIL5.g.vcf.gz \
 --variant BIL53.g.vcf.gz \
 --variant BIL6.g.vcf.gz \
 --variant BIL62.g.vcf.gz \
 --variant BIL7.g.vcf.gz \
 --variant BIL8.g.vcf.gz \
 --variant BIL9.g.vcf.gz \
 --variant BOR1.g.vcf.gz \
 --variant BOR11.g.vcf.gz \
 --variant BOR14.g.vcf.gz \
 --variant BOR15.g.vcf.gz \
 --variant BOR16.g.vcf.gz \
 --variant BOR18.g.vcf.gz \
 --variant BOR2.g.vcf.gz \
 --variant BOR20.g.vcf.gz \
 --variant BOR22.g.vcf.gz \
 --variant BOR23-Q.g.vcf.gz \
 --variant BOR24.g.vcf.gz \
 --variant BOR5.g.vcf.gz \
 --variant BOR9.g.vcf.gz \
 --variant BRO3-DUP.g.vcf.gz \
 --variant BRO3.g.vcf.gz \
 --variant BRO5.g.vcf.gz \
 --variant BRO6.g.vcf.gz \
 --variant BRO9.g.vcf.gz \
 --variant BUK15.g.vcf.gz \
 --variant BUK17.g.vcf.gz \
 --variant BUK5.g.vcf.gz \
 --variant BUK7.g.vcf.gz \
 --variant DEB1.g.vcf.gz \
 --variant DEB16.g.vcf.gz \
 --variant DEB2.g.vcf.gz \
 --variant DEB9.g.vcf.gz \
 --variant EFI1-2.g.vcf.gz \
 --variant EFI1-4.g.vcf.gz \
 --variant EFI1-7-DUP.g.vcf.gz \
 --variant EFI1-7.g.vcf.gz \
 --variant EFI1-8.g.vcf.gz \
 --variant EFI1-9.g.vcf.gz \
 --variant EFI2-10.g.vcf.gz \
 --variant EFI2-11.g.vcf.gz \
 --variant EFI2-12.g.vcf.gz \
 --variant EFI2-4.g.vcf.gz \
 --variant EFI2-5.g.vcf.gz \
 --variant EFI2-6.g.vcf.gz \
 --variant EFI2-8.g.vcf.gz \
 --variant EFI2-9.g.vcf.gz \
 --variant FRE3.g.vcf.gz \
 --variant FRE4.g.vcf.gz \
 --variant FRE5.g.vcf.gz \
 --variant FRE6.g.vcf.gz \
 --variant GOS1.g.vcf.gz \
 --variant GOS10.g.vcf.gz \
 --variant GOS12.g.vcf.gz \
 --variant GOS13.g.vcf.gz \
 --variant GOS14.g.vcf.gz \
 --variant GOS15.g.vcf.gz \
 --variant GOS16.g.vcf.gz \
 --variant GOS17.g.vcf.gz \
 --variant GOS18-Q.g.vcf.gz \
 --variant GOS3.g.vcf.gz \
 --variant GOS4.g.vcf.gz \
 --variant GOS7.g.vcf.gz \
 --variant GOS8.g.vcf.gz \
 --variant GOS9-Q.g.vcf.gz \
 --variant KUK2.g.vcf.gz \
 --variant KUK9.g.vcf.gz \
 --variant LAM10.g.vcf.gz \
 --variant LAM11.g.vcf.gz \
 --variant LAM3.g.vcf.gz \
 --variant LAM7.g.vcf.gz \
 --variant LAM8.g.vcf.gz \
 --variant LAN1-DUP.g.vcf.gz \
 --variant LAN1.g.vcf.gz \
 --variant LAN11.g.vcf.gz \
 --variant LAN16.g.vcf.gz \
 --variant LAN17.g.vcf.gz \
 --variant LAN19.g.vcf.gz \
 --variant LAN2.g.vcf.gz \
 --variant LAN25-H.g.vcf.gz \
 --variant LAN27.g.vcf.gz \
 --variant LAN3.g.vcf.gz \
 --variant LAN6.g.vcf.gz \
 --variant LAN7.g.vcf.gz \
 --variant LAN8.g.vcf.gz \
 --variant LAN9.g.vcf.gz \
 --variant LUB1-DUP.g.vcf.gz \
 --variant LUB1.g.vcf.gz \
 --variant LUB17.g.vcf.gz \
 --variant LUB19.g.vcf.gz \
 --variant LUB2.g.vcf.gz \
 --variant LUB21.g.vcf.gz \
 --variant LUB22.g.vcf.gz \
 --variant LUB31.g.vcf.gz \
 --variant LUB35.g.vcf.gz \
 --variant LUB38.g.vcf.gz \
 --variant LUB5.g.vcf.gz \
 --variant LUB6.g.vcf.gz \
 --variant LUB62-H.g.vcf.gz \
 --variant LUB7.g.vcf.gz \
 --variant LUB8.g.vcf.gz \
 --variant MEL17.g.vcf.gz \
 --variant MEL18.g.vcf.gz \
 --variant MEL19.g.vcf.gz \
 --variant MEL20.g.vcf.gz \
 --variant MEL22.g.vcf.gz \
 --variant MEL23.g.vcf.gz \
 --variant MEL24.g.vcf.gz \
 --variant MEL25.g.vcf.gz \
 --variant MEL26.g.vcf.gz \
 --variant MEL27.g.vcf.gz \
 --variant MEL28.g.vcf.gz \
 --variant MEL29.g.vcf.gz \
 --variant MEL30.g.vcf.gz \
 --variant MEL32.g.vcf.gz \
 --variant SIL1-4-DUP.g.vcf.gz \
 --variant SIL1-4.g.vcf.gz \
 --variant SIL1-5.g.vcf.gz \
 --variant SIL1-6.g.vcf.gz \
 --variant SIL1-8.g.vcf.gz \
 --variant SIL2-2.g.vcf.gz \
 --variant SIL2-6.g.vcf.gz \
 --variant SIL3-4.g.vcf.gz \
 --variant SIL3-5.g.vcf.gz \
 --variant SIL3-6.g.vcf.gz \
 --variant SIL4-2.g.vcf.gz \
 --variant SIL4-3.g.vcf.gz \
 --variant SIL4-6.g.vcf.gz \
 --variant SIL5-2.g.vcf.gz \
 --variant SIL5-7.g.vcf.gz \
 --variant STE1.g.vcf.gz \
 --variant STE12.g.vcf.gz \
 --variant STE2.g.vcf.gz \
 --variant STE20.g.vcf.gz \
 --variant STE22.g.vcf.gz \
 --variant STE24.g.vcf.gz \
 --variant STE25.g.vcf.gz \
 --variant STE3-H.g.vcf.gz \
 --variant STE31.g.vcf.gz \
 --variant STE33.g.vcf.gz \
 --variant STE4-DUP.g.vcf.gz \
 --variant STE4.g.vcf.gz \
 --variant STE5.g.vcf.gz \
 --variant STE8.g.vcf.gz \
 --variant STE9.g.vcf.gz \
 --variant STJ1-DUP.g.vcf.gz \
 --variant STJ1.g.vcf.gz \
 --variant STJ10.g.vcf.gz \
 --variant STJ11.g.vcf.gz \
 --variant STJ12.g.vcf.gz \
 --variant STJ14.g.vcf.gz \
 --variant STJ2.g.vcf.gz \
 --variant STJ24.g.vcf.gz \
 --variant STJ3.g.vcf.gz \
 --variant STJ4.g.vcf.gz \
 --variant STJ5.g.vcf.gz \
 --variant STJ6.g.vcf.gz \
 --variant STJ7.g.vcf.gz \
 --variant STJ8.g.vcf.gz \
 --variant STJ9.g.vcf.gz \
 --variant SVA1.g.vcf.gz \
 --variant SVA10.g.vcf.gz \
 --variant SVA11.g.vcf.gz \
 --variant SVA15.g.vcf.gz \
 --variant SVA16.g.vcf.gz \
 --variant SVA18.g.vcf.gz \
 --variant SVA3.g.vcf.gz \
 --variant SVA4.g.vcf.gz \
 --variant SVA5.g.vcf.gz \
 --variant SVA6.g.vcf.gz \
 --variant SVA7.g.vcf.gz \
 --variant SVA8-DUP.g.vcf.gz \
 --variant SVA8.g.vcf.gz \
 --variant SVA9.g.vcf.gz \
 --variant TON10.g.vcf.gz \
 --variant TON11.g.vcf.gz \
 --variant TON18.g.vcf.gz \
 --variant TON19.g.vcf.gz \
 --variant TON2.g.vcf.gz \
 --variant TON22-H.g.vcf.gz \
 --variant TON3.g.vcf.gz \
 --variant TON4.g.vcf.gz \
 --variant TON5.g.vcf.gz \
 --variant TON6.g.vcf.gz \
 --variant TON7.g.vcf.gz \
 --variant TON8.g.vcf.gz \
 --variant TON9.g.vcf.gz \
 --variant TRE1.g.vcf.gz \
 --variant TRE10.g.vcf.gz \
 --variant TRE11.g.vcf.gz \
 --variant TRE12.g.vcf.gz \
 --variant TRE15-Q.g.vcf.gz \
 --variant TRE16.g.vcf.gz \
 --variant TRE3.g.vcf.gz \
 --variant TRE4.g.vcf.gz \
 --variant TRE5.g.vcf.gz \
 --variant TRE6.g.vcf.gz \
 --variant TRE7.g.vcf.gz \
 --variant TRE8.g.vcf.gz \
 --variant TRE9.g.vcf.gz \
 -O Bark_beetle.g.vcf.gz

#13.3. GenotypeGVCFs - run over whole genome#
gatk --java-options "-Xmx40g" GenotypeGVCFs \
   -R SortedAndRenamed.clean.fasta \
   -V Bark_beetle.g.vcf.gz \
   -O Bark_beetle.vcf.gz

#13.4. GenotypeGVCFs was slow thus it had to be parallelized - by dividing file by contigs to separate files and then running it in paralell contig by contig#
#To run it in paralell - List of genomic intervals (.bed file) is needed#
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' SortedAndRenamed.clean.fasta.fai > SortedAndRenamed.clean.fasta.bed

#This list is later spread to separate files - one per contig#
cat SortedAndRenamed.clean.fasta.bed | while read line
do
echo $line > $(cut -d$'\t' -f1 <<< $line).bed
done

#Next this files were used to prepare list of them - over this list parallel version of Genotype GVCFs was running#
for file in IpsContig*.bed
	do
	echo $file >> Bed_list.txt
	done

#GenotypeGVCFs in parallel mode was iterating over list of .g.vcf files created in a previous step per contig#
for file in IpsContig*.g.vcf.gz
	do
	echo $file >> GVCF_list.txt
	done

#Merging single Contig VCFs into one big file - To be run after GenotypeGVCFs#
for file in IpsContig*.g.vcf.gz
	do
	genotyped=$(echo $file |cut -d"." -f 1,3,4)
	echo $genotyped >> VCFs.txt
	done

picard MergeVcfs \
	I=VCFs.txt \
	O=Bark_beetle_meged.vcf.gz

#14.FILTERING VARIANTS#
#Hard filtering steps#
gatk --java-options "-Xmx80g" VariantFiltration \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged.vcf \
	-O Bark_beetle_meged_Tech_filtering.vcf.gz \
	--filter-expression "FS > 60.0" \
	--filter-name "FS_filter" \
	--filter-expression "QD < 2.0" \
	--filter-name "QD_filter" \
	--filter-expression "MQ < 30.0" \
	--filter-name "MQ_filter" \
	--filter-expression "SOR > 3.0" \
	--filter-name "SOR_filter" \
	--filter-expression "MQRankSum < -12.5" \
	--filter-name "MQRankSum_filter" \
	--filter-expression "ReadPosRankSum < -8.0" \
	--filter-name "ReadPosRankSum_filter"

gatk --java-options "-Xmx80g" VariantFiltration \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_Tech_filtering.vcf.gz \
	-O Bark_beetle_meged_Tech_filtering_GT_filters.vcf.gz \
	-G-filter "GQ < 20.0" \
	-G-filter-name "GQlow" \
	-G-filter "DP < 8.0" \
	-G-filter-name "DPlow" \
	--set-filtered-genotype-to-no-call

gatk --java-options "-Xmx80g" VariantFiltration \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_Tech_filtering_GT_filters.vcf.gz \
	-O Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters.vcf.gz \
	--filter-expression "AF < 0.01" \
	--filter-name "MAF_L" \
	--filter-expression "AF > 0.99" \
	--filter-name "MAF_H" \
	--filter-expression "AN < 253.0" \
	--filter-name "MD"

#Remove 5 SNPs around indels SNPs bcftools#
bcftools filter -g 5 -O v -o Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels.vcf Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters.vcf.gz

#Now remove sites which were filtered - field FILTER different than PASS, OR are not biallelic SNPs#
#The ouput file from this step will serve as database of known SNP's#
gatk SelectVariants \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels.vcf \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	--exclude-filtered \
	--exclude-non-variants \
	--remove-unused-alternates \
	-O Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels_PASSING.vcf.gz

#15. BaseRecalibrator was run to detect systematic errors in base quality scores#
#BQSR - should be run in per lane mode#
#15.1. Since some filelds (especially platform and lane ID ie. -RGPL, -RGPU) were missing from initial .bam files BUT were reqired by GATK BQSR we introduced them with picard#
#All missing fileds were were identified from the read name/title in Fastq files#
#If one has one library for each sample running on one lane of a sequencing machine then you can make SM=LB=RGID=PU - from https://angus.readthedocs.io/en/2017/Read_group_info.html#
#Thus all necessary info/fields were added as in an exmple below#

for file in *_coordsorted_rmdups.bam
	do
	nameis=$(echo $file | cut -d'_' -f 1)
	picard AddOrReplaceReadGroups -INPUT $file -OUTPUT $nameis"_coordsorted_rmdups_addPLPU.bam" -RGID $nameis -RGLB $nameis -RGPL ILLUMINA -RGPU "HHWGFDSXY.4."$nameis -RGSM $nameis
	done

files=*_coordsorted_rmdups_addPLPU.bam
for plik in $files
  do
  nameis=$(echo $plik | cut -d'_' -f 1)
  echo $nameis
  samtools index $nameis"_coordsorted_rmdups_addPLPU.bam"
  done

#15.2. Base reclibrator, run separetely for indivudals from separte lanes#
gatk BaseRecalibrator \
    -R SortedAndRenamed.clean.fasta \
    -I AAS10_coordsorted_rmdups_addPLPU.bam \
    -I AAS18_coordsorted_rmdups_addPLPU.bam \
    -I AAS2_coordsorted_rmdups_addPLPU.bam \
    -I AAS21_coordsorted_rmdups_addPLPU.bam \
    -I AAS22_coordsorted_rmdups_addPLPU.bam \
    -I AAS24_coordsorted_rmdups_addPLPU.bam \
    -I AAS25-DUP_coordsorted_rmdups_addPLPU.bam \
    -I AAS25_coordsorted_rmdups_addPLPU.bam \
    -I AAS28_coordsorted_rmdups_addPLPU.bam \
    -I AAS29_coordsorted_rmdups_addPLPU.bam \
    -I AAS3_coordsorted_rmdups_addPLPU.bam \
    -I AAS30_coordsorted_rmdups_addPLPU.bam \
    -I AAS41_coordsorted_rmdups_addPLPU.bam \
    -I AAS6_coordsorted_rmdups_addPLPU.bam \
    -I ASA10_coordsorted_rmdups_addPLPU.bam \
    -I ASA11_coordsorted_rmdups_addPLPU.bam \
    -I ASA12_coordsorted_rmdups_addPLPU.bam \
    -I ASA13_coordsorted_rmdups_addPLPU.bam \
    -I ASA14_coordsorted_rmdups_addPLPU.bam \
    -I ASA2_coordsorted_rmdups_addPLPU.bam \
    -I ASA20_coordsorted_rmdups_addPLPU.bam \
    -I ASA21_coordsorted_rmdups_addPLPU.bam \
    -I ASA46_coordsorted_rmdups_addPLPU.bam \
    -I ASA47_coordsorted_rmdups_addPLPU.bam \
    -I ASA8_coordsorted_rmdups_addPLPU.bam \
    -I BAW148_coordsorted_rmdups_addPLPU.bam \
    -I BAW150_coordsorted_rmdups_addPLPU.bam \
    -I BAW152_coordsorted_rmdups_addPLPU.bam \
    -I BAW153_coordsorted_rmdups_addPLPU.bam \
    -I BAW154_coordsorted_rmdups_addPLPU.bam \
    -I BAW156_coordsorted_rmdups_addPLPU.bam \
    -I BAW158_coordsorted_rmdups_addPLPU.bam \
    -I BAW160_coordsorted_rmdups_addPLPU.bam \
    -I BAW162_coordsorted_rmdups_addPLPU.bam \
    -I BAW164_coordsorted_rmdups_addPLPU.bam \
    -I BAW165_coordsorted_rmdups_addPLPU.bam \
    -I BAW175_coordsorted_rmdups_addPLPU.bam \
    -I BAW178_coordsorted_rmdups_addPLPU.bam \
    -I EFI1-7-DUP_coordsorted_rmdups_addPLPU.bam \
    -I LAN1-DUP_coordsorted_rmdups_addPLPU.bam \
    -I LAN1_coordsorted_rmdups_addPLPU.bam \
    -I LAN11_coordsorted_rmdups_addPLPU.bam \
    -I LAN16_coordsorted_rmdups_addPLPU.bam \
    -I LAN17_coordsorted_rmdups_addPLPU.bam \
    -I LAN19_coordsorted_rmdups_addPLPU.bam \
    -I LAN2_coordsorted_rmdups_addPLPU.bam \
    -I LAN25-H_coordsorted_rmdups_addPLPU.bam \
    -I LAN27_coordsorted_rmdups_addPLPU.bam \
    -I LAN3_coordsorted_rmdups_addPLPU.bam \
    -I LAN6_coordsorted_rmdups_addPLPU.bam \
    -I LAN7_coordsorted_rmdups_addPLPU.bam \
    -I LAN8_coordsorted_rmdups_addPLPU.bam \
    -I LAN9_coordsorted_rmdups_addPLPU.bam \
    -I LUB62-H_coordsorted_rmdups_addPLPU.bam \
    -I SIL1-4-DUP_coordsorted_rmdups_addPLPU.bam \
    -I SIL1-4_coordsorted_rmdups_addPLPU.bam \
    -I SIL1-5_coordsorted_rmdups_addPLPU.bam \
    -I SIL1-6_coordsorted_rmdups_addPLPU.bam \
    -I SIL1-8_coordsorted_rmdups_addPLPU.bam \
    -I SIL2-2_coordsorted_rmdups_addPLPU.bam \
    -I SIL2-6_coordsorted_rmdups_addPLPU.bam \
    -I SIL3-4_coordsorted_rmdups_addPLPU.bam \
    -I SIL3-5_coordsorted_rmdups_addPLPU.bam \
    -I SIL3-6_coordsorted_rmdups_addPLPU.bam \
    -I SIL4-2_coordsorted_rmdups_addPLPU.bam \
    -I SIL4-3_coordsorted_rmdups_addPLPU.bam \
    -I SIL4-6_coordsorted_rmdups_addPLPU.bam \
    -I SIL5-2_coordsorted_rmdups_addPLPU.bam \
    -I SIL5-7_coordsorted_rmdups_addPLPU.bam \
    -I STE3-H_coordsorted_rmdups_addPLPU.bam \
    -I SVA1_coordsorted_rmdups_addPLPU.bam \
    -I SVA10_coordsorted_rmdups_addPLPU.bam \
    -I SVA11_coordsorted_rmdups_addPLPU.bam \
    -I SVA15_coordsorted_rmdups_addPLPU.bam \
    -I SVA16_coordsorted_rmdups_addPLPU.bam \
    -I SVA18_coordsorted_rmdups_addPLPU.bam \
    -I SVA3_coordsorted_rmdups_addPLPU.bam \
    -I SVA4_coordsorted_rmdups_addPLPU.bam \
    -I SVA5_coordsorted_rmdups_addPLPU.bam \
    -I SVA6_coordsorted_rmdups_addPLPU.bam \
    -I SVA7_coordsorted_rmdups_addPLPU.bam \
    -I SVA8-DUP_coordsorted_rmdups_addPLPU.bam \
    -I SVA8_coordsorted_rmdups_addPLPU.bam \
    -I SVA9_coordsorted_rmdups_addPLPU.bam \
    -I TON22-H_coordsorted_rmdups_addPLPU.bam \
    --known-sites Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels_PASSING.vcf.gz \
    -O lane1_recal_data_PU.table

gatk BaseRecalibrator \
    -R SortedAndRenamed.clean.fasta \
    -I ASA124_coordsorted_rmdups_addPLPU.bam \
    -I ASA125_coordsorted_rmdups_addPLPU.bam \
    -I ASA126_coordsorted_rmdups_addPLPU.bam \
    -I BAD1_coordsorted_rmdups_addPLPU.bam \
    -I BAD16_coordsorted_rmdups_addPLPU.bam \
    -I BAD17_coordsorted_rmdups_addPLPU.bam \
    -I BAD18_coordsorted_rmdups_addPLPU.bam \
    -I BAD3_coordsorted_rmdups_addPLPU.bam \
    -I BIL1-DUP_coordsorted_rmdups_addPLPU.bam \
    -I BOR23-Q_coordsorted_rmdups_addPLPU.bam \
    -I BRO3-DUP_coordsorted_rmdups_addPLPU.bam \
    -I BRO3_coordsorted_rmdups_addPLPU.bam \
    -I BRO5_coordsorted_rmdups_addPLPU.bam \
    -I BRO6_coordsorted_rmdups_addPLPU.bam \
    -I BUK15_coordsorted_rmdups_addPLPU.bam \
    -I BUK17_coordsorted_rmdups_addPLPU.bam \
    -I DEB16_coordsorted_rmdups_addPLPU.bam \
    -I DEB9_coordsorted_rmdups_addPLPU.bam \
    -I EFI1-2_coordsorted_rmdups_addPLPU.bam \
    -I EFI1-4_coordsorted_rmdups_addPLPU.bam \
    -I EFI1-7_coordsorted_rmdups_addPLPU.bam \
    -I EFI1-8_coordsorted_rmdups_addPLPU.bam \
    -I EFI1-9_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-10_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-11_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-12_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-4_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-5_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-6_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-8_coordsorted_rmdups_addPLPU.bam \
    -I EFI2-9_coordsorted_rmdups_addPLPU.bam \
    -I FRE3_coordsorted_rmdups_addPLPU.bam \
    -I FRE4_coordsorted_rmdups_addPLPU.bam \
    -I FRE5_coordsorted_rmdups_addPLPU.bam \
    -I FRE6_coordsorted_rmdups_addPLPU.bam \
    -I GOS1_coordsorted_rmdups_addPLPU.bam \
    -I GOS10_coordsorted_rmdups_addPLPU.bam \
    -I GOS12_coordsorted_rmdups_addPLPU.bam \
    -I GOS13_coordsorted_rmdups_addPLPU.bam \
    -I GOS14_coordsorted_rmdups_addPLPU.bam \
    -I GOS15_coordsorted_rmdups_addPLPU.bam \
    -I GOS16_coordsorted_rmdups_addPLPU.bam \
    -I GOS17_coordsorted_rmdups_addPLPU.bam \
    -I GOS18-Q_coordsorted_rmdups_addPLPU.bam \
    -I GOS3_coordsorted_rmdups_addPLPU.bam \
    -I GOS4_coordsorted_rmdups_addPLPU.bam \
    -I GOS7_coordsorted_rmdups_addPLPU.bam \
    -I GOS8_coordsorted_rmdups_addPLPU.bam \
    -I GOS9-Q_coordsorted_rmdups_addPLPU.bam \
    -I KUK9_coordsorted_rmdups_addPLPU.bam \
    -I LAM10_coordsorted_rmdups_addPLPU.bam \
    -I LAM11_coordsorted_rmdups_addPLPU.bam \
    -I LAM3_coordsorted_rmdups_addPLPU.bam \
    -I LAM7_coordsorted_rmdups_addPLPU.bam \
    -I LAM8_coordsorted_rmdups_addPLPU.bam \
    -I LUB1-DUP_coordsorted_rmdups_addPLPU.bam \
    -I LUB1_coordsorted_rmdups_addPLPU.bam \
    -I LUB17_coordsorted_rmdups_addPLPU.bam \
    -I LUB19_coordsorted_rmdups_addPLPU.bam \
    -I LUB2_coordsorted_rmdups_addPLPU.bam \
    -I LUB21_coordsorted_rmdups_addPLPU.bam \
    -I LUB22_coordsorted_rmdups_addPLPU.bam \
    -I LUB31_coordsorted_rmdups_addPLPU.bam \
    -I LUB35_coordsorted_rmdups_addPLPU.bam \
    -I LUB38_coordsorted_rmdups_addPLPU.bam \
    -I LUB5_coordsorted_rmdups_addPLPU.bam \
    -I LUB6_coordsorted_rmdups_addPLPU.bam \
    -I LUB7_coordsorted_rmdups_addPLPU.bam \
    -I LUB8_coordsorted_rmdups_addPLPU.bam \
    -I STE1_coordsorted_rmdups_addPLPU.bam \
    -I STE12_coordsorted_rmdups_addPLPU.bam \
    -I STE2_coordsorted_rmdups_addPLPU.bam \
    -I STE20_coordsorted_rmdups_addPLPU.bam \
    -I STE22_coordsorted_rmdups_addPLPU.bam \
    -I STE24_coordsorted_rmdups_addPLPU.bam \
    -I STE25_coordsorted_rmdups_addPLPU.bam \
    -I STE31_coordsorted_rmdups_addPLPU.bam \
    -I STE33_coordsorted_rmdups_addPLPU.bam \
    -I STE4-DUP_coordsorted_rmdups_addPLPU.bam \
    -I STE4_coordsorted_rmdups_addPLPU.bam \
    -I STE5_coordsorted_rmdups_addPLPU.bam \
    -I STE8_coordsorted_rmdups_addPLPU.bam \
    -I STE9_coordsorted_rmdups_addPLPU.bam \
    -I STJ1-DUP_coordsorted_rmdups_addPLPU.bam \
    -I TRE15-Q_coordsorted_rmdups_addPLPU.bam \
    --known-sites Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels_PASSING.vcf.gz \
    -O lane2_recal_data_PU.table

gatk BaseRecalibrator \
    -R SortedAndRenamed.clean.fasta \
    -I BIL10_coordsorted_rmdups_addPLPU.bam \
    -I BIL11_coordsorted_rmdups_addPLPU.bam \
    -I BIL20_coordsorted_rmdups_addPLPU.bam \
    -I BIL27_coordsorted_rmdups_addPLPU.bam \
    -I BIL34_coordsorted_rmdups_addPLPU.bam \
    -I BIL40_coordsorted_rmdups_addPLPU.bam \
    -I BIL5_coordsorted_rmdups_addPLPU.bam \
    -I BIL53_coordsorted_rmdups_addPLPU.bam \
    -I BIL6_coordsorted_rmdups_addPLPU.bam \
    -I BIL62_coordsorted_rmdups_addPLPU.bam \
    -I BIL7_coordsorted_rmdups_addPLPU.bam \
    -I BIL8_coordsorted_rmdups_addPLPU.bam \
    -I BIL9_coordsorted_rmdups_addPLPU.bam \
    -I BOR1_coordsorted_rmdups_addPLPU.bam \
    -I BOR11_coordsorted_rmdups_addPLPU.bam \
    -I BOR14_coordsorted_rmdups_addPLPU.bam \
    -I BOR15_coordsorted_rmdups_addPLPU.bam \
    -I BOR16_coordsorted_rmdups_addPLPU.bam \
    -I BOR18_coordsorted_rmdups_addPLPU.bam \
    -I BOR2_coordsorted_rmdups_addPLPU.bam \
    -I BOR20_coordsorted_rmdups_addPLPU.bam \
    -I BOR22_coordsorted_rmdups_addPLPU.bam \
    -I BOR24_coordsorted_rmdups_addPLPU.bam \
    -I BOR5_coordsorted_rmdups_addPLPU.bam \
    -I BOR9_coordsorted_rmdups_addPLPU.bam \
    -I BRO9_coordsorted_rmdups_addPLPU.bam \
    -I BUK5_coordsorted_rmdups_addPLPU.bam \
    -I BUK7_coordsorted_rmdups_addPLPU.bam \
    -I DEB1_coordsorted_rmdups_addPLPU.bam \
    -I DEB2_coordsorted_rmdups_addPLPU.bam \
    -I KUK2_coordsorted_rmdups_addPLPU.bam \
    -I MEL17_coordsorted_rmdups_addPLPU.bam \
    -I MEL18_coordsorted_rmdups_addPLPU.bam \
    -I MEL19_coordsorted_rmdups_addPLPU.bam \
    -I MEL20_coordsorted_rmdups_addPLPU.bam \
    -I MEL22_coordsorted_rmdups_addPLPU.bam \
    -I MEL23_coordsorted_rmdups_addPLPU.bam \
    -I MEL24_coordsorted_rmdups_addPLPU.bam \
    -I MEL25_coordsorted_rmdups_addPLPU.bam \
    -I MEL26_coordsorted_rmdups_addPLPU.bam \
    -I MEL27_coordsorted_rmdups_addPLPU.bam \
    -I MEL28_coordsorted_rmdups_addPLPU.bam \
    -I MEL29_coordsorted_rmdups_addPLPU.bam \
    -I MEL30_coordsorted_rmdups_addPLPU.bam \
    -I MEL32_coordsorted_rmdups_addPLPU.bam \
    -I STJ1_coordsorted_rmdups_addPLPU.bam \
    -I STJ10_coordsorted_rmdups_addPLPU.bam \
    -I STJ11_coordsorted_rmdups_addPLPU.bam \
    -I STJ12_coordsorted_rmdups_addPLPU.bam \
    -I STJ14_coordsorted_rmdups_addPLPU.bam \
    -I STJ2_coordsorted_rmdups_addPLPU.bam \
    -I STJ24_coordsorted_rmdups_addPLPU.bam \
    -I STJ3_coordsorted_rmdups_addPLPU.bam \
    -I STJ4_coordsorted_rmdups_addPLPU.bam \
    -I STJ5_coordsorted_rmdups_addPLPU.bam \
    -I STJ6_coordsorted_rmdups_addPLPU.bam \
    -I STJ7_coordsorted_rmdups_addPLPU.bam \
    -I STJ8_coordsorted_rmdups_addPLPU.bam \
    -I STJ9_coordsorted_rmdups_addPLPU.bam \
    -I TON10_coordsorted_rmdups_addPLPU.bam \
    -I TON11_coordsorted_rmdups_addPLPU.bam \
    -I TON18_coordsorted_rmdups_addPLPU.bam \
    -I TON19_coordsorted_rmdups_addPLPU.bam \
    -I TON2_coordsorted_rmdups_addPLPU.bam \
    -I TON3_coordsorted_rmdups_addPLPU.bam \
    -I TON4_coordsorted_rmdups_addPLPU.bam \
    -I TON5_coordsorted_rmdups_addPLPU.bam \
    -I TON6_coordsorted_rmdups_addPLPU.bam \
    -I TON7_coordsorted_rmdups_addPLPU.bam \
    -I TON8_coordsorted_rmdups_addPLPU.bam \
    -I TON9_coordsorted_rmdups_addPLPU.bam \
    -I TRE1_coordsorted_rmdups_addPLPU.bam \
    -I TRE10_coordsorted_rmdups_addPLPU.bam \
    -I TRE11_coordsorted_rmdups_addPLPU.bam \
    -I TRE12_coordsorted_rmdups_addPLPU.bam \
    -I TRE16_coordsorted_rmdups_addPLPU.bam \
    -I TRE3_coordsorted_rmdups_addPLPU.bam \
    -I TRE4_coordsorted_rmdups_addPLPU.bam \
    -I TRE5_coordsorted_rmdups_addPLPU.bam \
    -I TRE6_coordsorted_rmdups_addPLPU.bam \
    -I TRE7_coordsorted_rmdups_addPLPU.bam \
    -I TRE8_coordsorted_rmdups_addPLPU.bam \
    -I TRE9_coordsorted_rmdups_addPLPU.bam \
    --known-sites Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels_PASSING.vcf.gz \
    -O lane3_recal_data_PU.table

#15.3. Apply BQSR - run in separate folders for lanes - requires coping bams to these folders first - IMPORTANT: give proper recal table to each folder/ set of individuals!#
for file in *_coordsorted_rmdups_addPLPU.bam
	do
	echo $file
	nameis=$(echo $file | cut -d'_' -f 1)
	echo $nameis
	gatk ApplyBQSR \
	-R SortedAndRenamed.clean.fasta \
	-I $file \
	-O $nameis"_recalibrated.bam" \
	--bqsr-recal-file lane3_recal_data_PU.table
	done

#16. !!! After BQSR was applied, RUN: HaplotypeCaller (on recalibrated bams!), CombineGVCFs, GenotypeGVCFs and Picard MergeVcfs!#

#17. Variant filtering#
#17.1. Perform VQSR - Variant Quality Score Recalibration#
gatk VariantRecalibrator \
   -R SortedAndRenamed.clean.fasta \
   -V Bark_beetle_meged.vcf.gz \
   --resource:dbsnp,known=true,training=true,truth=true Bark_beetle_meged_Tech_filtering_GT_filters_VARIANT_filters_no_sites_around_indels_PASSING.vcf.gz \
   -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
   -mode SNP \
   -O output.recal \
   --tranches-file output.tranches

gatk ApplyVQSR \
   -R SortedAndRenamed.clean.fasta \
   -V Bark_beetle_meged.vcf.gz \
   -O Bark_beetle_meged_VQSR.vcf.gz \
   --truth-sensitivity-filter-level 99.0 \
   --tranches-file output.tranches \
   --recal-file output.recal \
   -mode SNP

#17.2 Calculate how namy SNPs were filtered by VQSR filtering step#
vcftools --gzvcf Bark_beetle_meged.vcf.gz --FILTER-summary --out Raw
vcftools --gzvcf Bark_beetle_meged_VQSR.vcf.gz --FILTER-summary --out VQSR_filtering

#17.3. Remove 5 SNPs around indels SNPs bcftools#
bcftools filter -g 5 -O z -o Bark_beetle_meged_VQSR_no_sites_around_indels.vcf.gz Bark_beetle_meged_VQSR.vcf.gz
vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels.vcf.gz --FILTER-summary --out VQSR_filtering_no_sites_around_indels

gatk IndexFeatureFile \
     -I Bark_beetle_meged_VQSR_no_sites_around_indels.vcf.gz

#17.4. Remove sites which were filtered or are not biallelic SNPs#
gatk SelectVariants \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels.vcf.gz \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	--exclude-filtered \
	--exclude-non-variants \
	--remove-unused-alternates \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs.vcf.gz

vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs.vcf.gz --FILTER-summary --out VQSR_just_Biallelic_SNPs

#17.5. Hard filtering with filtering variants with excessive overall coverage (mean + 1 SD = 9319), and variants that exhibiting significant heterozygote excess (ExcessHet > 54.69)#
gatk --java-options "-Xmx80g" VariantFiltration \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs.vcf.gz \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering.vcf.gz \
	--filter-expression "FS > 60.0" \
	--filter-name "FS_filter" \
	--filter-expression "QD < 2.0" \
	--filter-name "QD_filter" \
	--filter-expression "MQ < 30.0" \
	--filter-name "MQ_filter" \
	--filter-expression "SOR > 3.0" \
	--filter-name "SOR_filter" \
	--filter-expression "DP > 9319" \
	--filter-name "DP_H" \
	--filter-expression "ExcessHet > 54.69" \
	--filter-name "ExcessHet"

vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering.vcf.gz --FILTER-summary --out SNPs_after_hard_filtering

#17.6. Filtering poor quality and low coverage genotypes#
gatk --java-options "-Xmx80g" VariantFiltration \
 -R SortedAndRenamed.clean.fasta \
 -V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering.vcf.gz \
 -O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters.vcf.gz \
 -G-filter "GQ < 20.0" \
 -G-filter-name "GQlow" \
 -G-filter "DP < 8.0" \
 -G-filter-name "DPlow" \
 --set-filtered-genotype-to-no-call

#17.7. Variants with too much (more than 50%) Missing data removal#
gatk --java-options "-Xmx80g" VariantFiltration \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters.vcf.gz \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD.vcf.gz \
	--filter-expression "AN < 253.0" \
	--filter-name "MD" \

vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD.vcf.gz --FILTER-summary --out SNPs_after_hard_filtering_and_MD_filtering

#17.8. Remove sites which were filtered on previous steps or are not biallelic SNPs.#
gatk SelectVariants \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD.vcf.gz \
	--select-type-to-include SNP \
	--restrict-alleles-to BIALLELIC \
	--exclude-filtered \
	--exclude-non-variants \
	--remove-unused-alternates \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING.vcf.gz

vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING.vcf.gz --FILTER-summary --out SNPs_after_hard_filtering_and_MD_filtering_PASSING
#17.9. Cut repetitive region#
bcftools view -O z -T ^Ref_coord.bed Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING.vcf.gz -o Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats.vcf.gz
gatk IndexFeatureFile \
     -I Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats.vcf.gz

vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats.vcf.gz --FILTER-summary --out SNPs_MD_filtered_after_repeat_masking

#18. Concordance test#
gatk SelectVariants \
 -R SortedAndRenamed.clean.fasta \
 -V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats.vcf.gz \
 -O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_COMP.vcf \
 -sn AAS25-DUP \
 -sn BRO3-DUP \
 -sn EFI1-7-DUP \
 -sn LAN1-DUP \
 -sn LUB1-DUP \
 -sn SIL1-4-DUP \
 -sn STE4-DUP \
 -sn STJ1-DUP \
 -sn SVA8-DUP

gatk SelectVariants \
 -R SortedAndRenamed.clean.fasta \
 -V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats.vcf.gz \
 -O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_EVAL.vcf \
 -sn AAS25 \
 -sn BRO3 \
 -sn EFI1-7 \
 -sn LAN1 \
 -sn LUB1 \
 -sn SIL1-4 \
 -sn STE4 \
 -sn STJ1 \
 -sn SVA8

#Before actual concordance test one need to manually change individual id's in "...COMP.vcf" - "DUP" had to be erased - so the names of individuals in COMP & EVAL will be exactly the same# 
nano Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_COMP.vcf

#Actual concordance asessment#
java -jar /home/piotr.zielinski/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
 -T GenotypeConcordance \
 -R SortedAndRenamed.clean.fasta \
 -eval Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_EVAL.vcf \
 -comp Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_COMP.vcf \
 -o Concordance_NO_MAF_filter.grp

#19. Removing duplicated sampels from the dataset#
gatk SelectVariants \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats.vcf.gz \
	-xl-sn AAS25-DUP \
	-xl-sn BRO3-DUP \
	-xl-sn EFI1-7-DUP \
	-xl-sn LAN1-DUP \
	-xl-sn LUB1-DUP \
	-xl-sn SIL1-4-DUP \
	-xl-sn STE4-DUP \
	-xl-sn STJ1-DUP \
	-xl-sn SVA8-DUP \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP.vcf.gz

#20. Initial whole genome PCA to asess if there is some structure in the data#
plink --vcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP.vcf.gz --pca header --allow-extra-chr

#21. Check if and which individuals are too similar to each other than expected by chance#
plink --vcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP.vcf.gz --genome --allow-extra-chr

#22. Removing sampels which were too similar to each other (as evidenced by IBS test in plink)#
gatk SelectVariants \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP.vcf.gz \
	-xl-sn STJ5 \
	-xl-sn BOR20 \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP_no_REL.vcf.gz

#23. Removing unused alternates - after excuding some individuals there could be some#
gatk SelectVariants \
	-R SortedAndRenamed.clean.fasta \
	-V Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP_no_REL.vcf.gz \
	--exclude-non-variants \
	--remove-unused-alternates \
	-O Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP_no_REL_PASSING.vcf.gz

#24. Extract single contig vcf and perform PCA contig by contig#
tabix --list-chroms Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP_no_REL_PASSING.vcf.gz > IpsContigs.txt

while IFS= read -r line; do
  vcftools --gzvcf Bark_beetle_meged_VQSR_no_sites_around_indels_Biallelic_SNPs_Hard_filtering_GT_filters_MD_PASSING_No_Repeats_no_DUP_no_REL_PASSING.vcf.gz --chr $line --recode --recode-INFO-all --out VCF_$line;
done < IpsContigs.txt

#Single contig PCA#
plink --vcf VCF_IpsContig1.recode.vcf --pca header --allow-extra-chr

#OR more automatic version#
for i in *recode.vcf
  do
  ind_1=$(echo $i |cut -d"." -f 1)
  ind_2=$(echo $ind_1 |cut -d"_" -f 2)
  plink --vcf $i --pca header --allow-extra-chr --out $ind_2
  done

#BAW152 & BAW154 the two most covered (coverage > mean*3) sampels look odd on the PCA plots - therefore they were also removed#
