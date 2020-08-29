#program to estimate number of rDNA copies

import sys
from scipy import stats
import numpy as np
from sklearn.linear_model import LinearRegression

#set up sliding window function
def slideFunct(data, window, step):
	total = len(data)
	spots = []
	for i in range (0, total-window, step):
		spots.append(i)
	result = []
	for i in range (len(spots)):
		result.append(np.mean(data[spots[i]:spots[i]+window]))
	return(result)

#depth file parser
def parser_dfile(depth_file):
	info = []
	with open(depth_file, 'r') as file_handle:                     #open reference fasta file for reading
		for line in file_handle:                        #loop of the lines in the config file
			col = 0
			char_list = ''
			for index in range(len(line)):                  #loop to parse each line
				if col < 2:				#skip two colums
					if (line[index].isspace()):     #if space or tab
						col += 1
					continue
				elif (line[index].isdigit()):		#if digit
					char_list += line[index]        #add value to the character list
				else:					#if not digit
					info.append(char_list)          #inserting value in info
	return info

##examine the read depth across the whole of Chromosome XII
def examine_whole_chromosome(whole_genome_depth, rdnastart, rdnaend, bufer, window, step):
	wgsdata12 = parser_dfile(whole_genome_depth)		#parsing the depth file
	wgsdata12 = np.array(wgsdata12).astype(np.int)
	startdata12 = []
	for j in range(rdnastart-bufer):
		startdata12.append(wgsdata12[j])
	enddata12 = []
	for j in range(rdnaend+bufer, len(wgsdata12)):
		enddata12.append(wgsdata12[j])
	notrdna12 = np.concatenate((startdata12, enddata12))
	sw12 = slideFunct(notrdna12, window, step)
	length = len(sw12)
	sw12index = []
	for j in range(length):
		sw12index.append(j)
	slope, intercept, r_value, p_value, std_err = stats.linregress(sw12index,sw12)
	return intercept

#examine the read depth across the rDNA single unit, forcing the regression line to be horizontal
def examine_rdna_unit(rDNA_depth, bufer, window, step):
	rdnadata = parser_dfile(rDNA_depth)
	rdnadata = np.array(rdnadata).astype(np.int)
	for j in range(bufer):
		rdnadata[len(rdnadata)-2*bufer+j] += rdnadata[j]              #adding right shoulder
		rdnadata[bufer+j] += rdnadata[len(rdnadata)-bufer+j]        #adding left shoulder
		rdnadata2 = []
		for j in range(bufer+1, len(rdnadata)-bufer):
			rdnadata2.append(rdnadata[j])
		sw12 = slideFunct(rdnadata2, window, step)
	return sw12

#main function
def examine(whole_genome_depth, rdnastart, rdnaend, buffer1, rDNA_depth, buffer2, window, step):
	#examine the read depth across the whole of Chromosome XII
	intercept = examine_whole_chromosome(whole_genome_depth, rdnastart, rdnaend, buffer1, window, step)
	#examine the read depth across the rDNA single unit, forcing the regression line to be horizontal
	sw12 = examine_rdna_unit(rDNA_depth, buffer2, window, step)
	return int(round(np.mean(sw12)/intercept))

