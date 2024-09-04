for file in Inversion_list_*.bed
	do
	echo $file >> File_list.txt
	python ORs_in_Inversions_2.py $file >> ORs_in_Inversions.txt
	done
