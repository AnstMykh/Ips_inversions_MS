#done per populations (bams.txt is a filelist containing the full path for each bam file with one filename per row. )
./angsd -bam bams.txt -anc SortedAndRenamed.clean.fasta -gl 1 -P 10 -out BAW_10 -doSaf 1

#
list=$(ls *.saf.idx)


for i in $list
    do
    out=$(echo $i | cut -d '.' -f 1)
    realSFS $i -fold 1 -P 60 > $out.sfs
    done
# If you don't have the ancestral state, you can instead estimate the folded SFS
# ...This is done by supplying the -anc with the reference genome and applying -fold 1 to realSFS.



# Calculate the thetas for each site. The output from the above command are two files out.thetas.gz and out.thetas.idx
for i in $list
    do
    out=$(echo $i | cut -d '.' -f 1)
    realSFS saf2theta $i -outname $out -sfs $out.sfs -fold 1
    done

# Estimate Tajimas D and other statistics
list=$(ls *.saf.idx| cut -d '.' -f 1-2)

for i in $list
    do
    thetaStat do_stat $i.thetas.idx -win 50000 -step 50000  -outnames $i.50kb.thetasWindow.gz 
    thetaStat do_stat $i.thetas.idx -win 100000 -step 100000  -outnames $i.100kb.thetasWindow.gz 
    thetaStat do_stat $i.thetas.idx -win 200000 -step 200000  -outnames $i.200kb.thetasWindow.gz 
    done
