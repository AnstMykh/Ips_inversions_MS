# Written by Piotr Zielinski #
# 21.04.2023 #
# Script for testing if there is an enrichment in ORs in inverted regions.#
# Prepare the genome parameters and parmaeters of the inversions ie. genome length and inversion lenghts#
# DESIGN & ASSUMPTION: ORs are located in always the same places but Inversion places can change#
# Inversions should be sorted from longest to shortest for convinience#

import random

def rand():
	counter = 1
	while counter <= 10000:
		out = open ('Inversion_list_'+str(counter)+'.bed', 'w')
		genome_bed = 'Genome2.bed'
		rem_regions = []
		inversions = open ('Inversions.bed', 'r')

		with open(genome_bed) as a:
			for line in a:
				if line != "":
					line = line.strip()
					genome_data = line.split('\t')
					genome_start = genome_data[1]
					genome_start_num = int(genome_start)
					genome_end = genome_data[2]
					genome_end_num = int(genome_end)
					rem_region_first = [genome_start_num, genome_end_num]
					rem_regions.append(rem_region_first)

		for line in inversions.readlines():
			line_strip = line.strip()
			fragment = line_strip.split('\t')
			inv_name = fragment[0]
			start = fragment[1]
			end = fragment[2]
			start_num = int(start)
			end_num = int(end)
			inv_len = end_num - start_num
			interval = random.sample(rem_regions, 1)
			interval_coords = interval[0]
			interval_start = interval_coords[0]
			interval_end = interval_coords[1]
			interval_len = interval_end - interval_start
			if inv_len < interval_len:
				max_inv_start_in_intrval = interval_end - inv_len
				inv_start = random.randrange(interval_start, max_inv_start_in_intrval)
				inv_end = inv_start + inv_len
				out.write(str(inv_name)+'\t'+str(inv_start)+'\t'+str(inv_end)+'\n')
				rem_region_c = [interval_start, inv_start]
				rem_regions.append(rem_region_c)
				rem_region_d = [inv_end, interval_end]
				rem_regions.append(rem_region_d)
				rem_regions.remove(interval_coords)
				print ("Inversion placed sucessfully")
			if inv_len > interval_len:
				print ("Interval length is too short to fit desired inversion - retrying")
				iterator = 0
				while iterator < 1:
					interval_2 = random.sample(rem_regions, 1)
					interval_2_coords = interval_2[0]
					interval_2_start = interval_2_coords[0]
					interval_2_end = interval_2_coords[1]
					interval_2_len = interval_2_end - interval_2_start
					if inv_len < interval_2_len:
						max_inv_start_in_intrval_2 = interval_2_end - inv_len
						inv_start = random.randrange(interval_2_start, max_inv_start_in_intrval_2)
						inv_end = inv_start + inv_len
						out.write(str(inv_name)+'\t'+str(inv_start)+'\t'+str(inv_end)+'\n')
						rem_region_e = [interval_2_start, inv_start]
						rem_regions.append(rem_region_e)
						rem_region_f = [inv_end, interval_2_end]
						rem_regions.append(rem_region_f)
						rem_regions.remove(interval_2_coords)
						print ("Inversion placed sucessfully")
						iterator += 1
					if inv_len > interval_2_len:
						print ("Interval length is too short to fit desired inversion - retrying")
						continue
		print ("Iteration "+str(counter)+" done")
		counter += 1


rand()
