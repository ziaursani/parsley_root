#comparing vcf files

import sys
from os import path

def parser(VCF, f):
	with open(VCF, 'r') as file_handle:        #open file for readin
		info = []			#complete information array
		locus = []                      #coordinates of all loci
		ref = []                        #reference allele array
		alt = []                        #alternate allele array
		if f == 's':			#if simulation file
                	cnv = []                        #copy numbers of this variant (array)
		else:				#if prediction file
               		AO = []                        #Number of alternate allele
                	DP = []                         #depth of reads
		for line in file_handle:
			space_done =  [False for i in range(8)]   #space not arrived yet
			locus_temp = ''                 #coordinate of one locus
			ref_temp = ''                   #one reference allele
			alt_temp = ''                   #one alternate allele
			if f == 's':			#if simulation file
                		cnv_temp = ''                   #copy numbers of one variant
			else:				#if prediction file
                        	AO_temp = ''                   #frequency of one variant
                        	DP_temp = ''
			if line.startswith('#'):            #if header
				continue
			for index in range(len(line)):
				if not space_done[0]:
					if line[index].isspace():
						space_done[0] = True		#first space arrived
				elif not space_done[1]:
					if line[index].isdigit():
						locus_temp += line[index]
					else:
						space_done[1] = True		#second space arrived
				elif not space_done[2]:
                                        if line[index].isspace():
                                                space_done[2] = True            #third space arrived
				elif not space_done[3]:
					if line[index].isalpha():
						ref_temp += line[index]
					else:
						space_done[3] = True		#fourth space arrived
				elif not space_done[4]:
                                        if line[index].isalpha():
                                                alt_temp += line[index]
                                        else:
                                                space_done[4] = True            #fifth space arrived
				elif f == 's':					#if simulation file
					if line[index].isdigit():
                                		cnv_temp += line[index]

					else:
						locus.append(locus_temp)
						ref.append(ref_temp)
						alt.append(alt_temp)
						cnv.append(cnv_temp)
						break
				elif not space_done[5]:
                                        if line[index] == '=' and line[index-1] == 'O' and line[index-2] == 'A' and line[index-3] == ';':
                                                space_done[5] = True
				elif not space_done[6]:
                                        if line[index].isdigit():
                                                AO_temp += line[index]
                                        else:
                                                space_done[6] = True
				elif not space_done[7]:
                                        if line[index] == '=' and line[index-1] == 'P' and line[index-2] == 'D' and line[index-3] == ';':
                                                space_done[7] = True
				elif line[index].isdigit():
                                        DP_temp += line[index]
				else:
                                        locus.append(locus_temp)
                                        ref.append(ref_temp)
                                        alt.append(alt_temp)
                                        AO.append(AO_temp)
                                        DP.append(DP_temp)
                                        break
		info.append(locus)
		info.append(ref)
		info.append(alt)
		if f == 's':				#if simulation file
			info.append(cnv)
		else:					#if prediction file
			info.append(AO)
			info.append(DP)
	return info

#comparing two vcfs
def VCF_comparer(sim_info, prd_info, cn):
	numvariants_sim = len(sim_info[0])
	numvariants_prd = len(prd_info[0])
	allele_done =  [False for i in range(numvariants_prd)]  
	comp_info = []					#comparative vcf info
	for i in range(numvariants_sim):		#variants from simulation vcf
		for j in range(3):
			comp_info.append(sim_info[j][i])
		comp_info.append(str(format(float(sim_info[3][i])/float(cn), '.3f')))	#frequency in simulated vcf
		match = False
		for j in range(numvariants_prd):		#matching variants from predictive vcf
			if prd_info[0][j] == sim_info[0][i] and  prd_info[1][j] == sim_info[1][i] and  prd_info[2][j] == sim_info[2][i]:
				comp_info.append(str(format(float(prd_info[3][j])/float(prd_info[4][j]), '.3f'))) #frequency in predictive vcf
				comp_info.append('Exact Match')
				allele_done[j] = True
				match = True
				break
			elif int(prd_info[0][j]) >= int(sim_info[0][i]) - 5 and int(prd_info[0][j]) <= int(sim_info[0][i]) and len(sim_info[1][i]) != len(sim_info[2][i]) and len(prd_info[1][j]) - len(prd_info[2][j]) == len(sim_info[1][i]) - len(sim_info[2][i]):
				comp_info.append(str(format(float(prd_info[3][j])/float(prd_info[4][j]), '.3f')))
				comp_info.append('Left Aligned')
				allele_done[j] = True
				match = True
				break
		if not match:
			comp_info.append('0.000')	#frequency in predictive vcf of non matching variants
			comp_info.append('Not Found')
	for i in range(numvariants_prd):		#variants from predictive vcf which did not match
		if allele_done[i]:			#if already matched
			continue
		for j in range(3):
                        comp_info.append(prd_info[j][i])
		comp_info.append('0.000')			#this did not match with the simulated vcf
		comp_info.append(str(format(float(prd_info[3][j])/float(prd_info[4][j]), '.3f'))) #frequency in predictive vcf
		comp_info.append('False Positive')
		
	return comp_info

#function to exit with error message
def exit_program(error_message):
        print(error_message)
        sys.exit(2)

#function to verify input
def input_check(args, i, Type):
        if(i+1 == len(args)):
                exit_program('Error1: Value missing: '+args[i])
        if(args[i+1].startswith('-')):
                exit_program('Error1: Value missing: '+args[i])
        if Type == 'int':
                if not args[i+1].isdigit():
                        exit_program('Error2: Not a number: '+args[i])
        elif Type == 'float':
                if not args[i+1].replace('.','',1).isdigit():
                        exit_program('Error3: Not numeric: '+args[i])
                if float(args[i+1]) < 0.0 or float(args[i+1]) > 1.0:
                        exit_program('Error4: Out of range: '+args[i])
        elif Type == 'str':
                if not path.exists(args[i+1]):
                        exit_program('Error5: Does not exist: '+args[i])
        return args[i+1], True	

#class to take command line arguments
class compare_variants():
	def __init__(self, args, args_len):
		global simulated_vcf, predictive_vcf, comparative_vcf, copy_number, python
		sim_vcf = False
		prd_vcf = False
		com_vcf = False
		cn = False
		for i in range(args_len-1):
			if args[i].startswith('-'):     #this is an argument attribute
				if args[i] == '--simulated_vcf' or args[i] == '-I':
                                        simulated_vcf, sim_vcf = input_check(args, i, 'str')                     #reference file
				elif args[i] == '--predictive_vcf' or args[i] == '-i':
                                        predictive_vcf, prd_vcf = input_check(args, i, 'str')     #value will only be assigned after input check
				elif args[i] == '--comparative_vcf' or args[i] == '-o':
					comparative_vcf = args[i+1]
					com_vcf = True
				elif args[i] == '--copy_number' or args[i] == '-n':
					copy_number, cn = input_check(args, i, 'int')
				else:
                                        exit_program('Error0: Argument not recognised '+args[i])
                #error and warning list
		if not sim_vcf:
                        exit_program('Error6: Missing argument: --simulated_vcf or -I.')
		if not prd_vcf:
                        exit_program('Error6: Missing argument: --predictive_vcf or -i.')
		if not com_vcf:
                        exit_program('Error6: Missing argument: --comparative_vcf or -o.')
		if not cn:
			exit_program('Error6: Missing argument: --copy_number or -n.')
		python = 'python'+str(sys.version_info.major)   #python version
	def execute(self):
                main()    

def main():
	sim_info = parser(simulated_vcf, 's')
	prd_info = parser(predictive_vcf, 'p')
	comp_info = VCF_comparer(sim_info, prd_info, copy_number)
	with open(comparative_vcf, 'w') as file_handle:
		file_handle.write('Locus'.ljust(25) + 'Ref'.ljust(25) + 'alt'.ljust(25)  + 'frq_sim'.ljust(25) + 'frq_prd'.ljust(25) + 'Comment' + '\n')
		for i in range(len(comp_info)):
			if (int((i+1)/6) < (i+1)/6):
				file_handle.write(comp_info[i].ljust(25))
			else:
				file_handle.write(comp_info[i] + '\n')

