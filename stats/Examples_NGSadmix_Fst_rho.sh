#NGSadmix

for rep in 1 2 3 4 5 6 7 8 9 10; do for i in 1 2 3 4 5; do NGSadmix -likes IpsContigNonInverted.BEAGLE.PL -K $i -P 10 -o IpsContigNonInverted_K$i"_rep"$rep -minMaf 0.05 -maxiter 10000; done; done


## example Fst in widnows; txt files with individuals names
vcftools --vcf data.vcf --weir-fst-pop NORTH.txt --weir-fst-pop SOUTH.txt --fst-window-size 100000 --fst-window-step 20000 --out NORTH.SOUTH

#rho calculation for example datafile
./complete
#settings n =40, low theta from angsd calculations, 4Ner =100, grid 101
./interval -loc SN.Contigs$contig.$i.ldhat.locs -seq SN.Contigs$contig.$i.ldhat.sites -lk SN_n40_complete_lk.txt -its 2000000 -bpen 5 -samp 5000
./stat -input rates.txt -burn 20 -loc SN.Contigs$contig.$i.ldhat.locs
#res.txt renamed SN.Contigs$contig.$i.ldhat.2MLN.pen5.samp5.COMPLETE.res.txt

#calculate rho in windows

perl windows.rho.pl res.txt res.txt.windowed.100kb 100

#or in sliding window
perl rho_windowed_slide.pl  res.txt.windowed.100kb.20kb 100 20
