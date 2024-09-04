# Written by Piotr Zielinski #
# 22.04.2023 #
# This script will check how many ORs is located within inverted regions#
# First prepare bed file with ORs - assumend to be always the same#
# Second prepare bed file with Inversion coordinates  - many such files are prepared automatically with Random_distribution_of_inversions.py #
# This script will take one command line argument - it should contain name of the bed file w ith inversion coordinates #
# Thus script can be run automatically with ORs_in_inversions_RUNNER.sh #

import sys

def rand():

	ORs_list = 'ORs.bed'
	ORs = set()

	with open(ORs_list) as p:
		for line in p:
			line = line.strip()
			if line != "":
				ORs.add(line)

	counter = 0
	#iterator = 0
	Inversions = []
	Inversions_bed = sys.argv[1]
	with open(Inversions_bed) as a:
		for line in a:
			if line != "":
				line = line.strip()
				inv_data = line.split('\t')
				inv_start = inv_data[1]
				inv_start_num = int(inv_start)
				inv_end = inv_data[2]
				inv_end_num = int(inv_end)
				for line in ORs:
					line_strip = line.strip()
					fragment = line_strip.split('\t')
					OR = fragment[0]
					start = fragment[1]
					end = fragment[2]
					start_num = int(start)
					end_num = int(end)
					if inv_start_num <= start_num <= inv_end_num or inv_start_num <= end_num <= inv_end_num or (inv_start_num <= start_num and end_num <= inv_end_num):
						counter += 1
						#print ("OR in inversion")
					#else:
						#iterator += 1
						#print ("OR outside inversion")
	print (counter)
	#print (iterator)
	
rand()
