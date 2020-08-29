import getopt, sys, os, subprocess, random
from os import path
from itertools import islice
from fasta_trimmer import trim_fasta
from copyno import examine

def max_line_length(infile):
	max_length = 0
	max_line_length = ''
	with open(infile, 'r') as inf:
		while True:
			fastq_line_block = list(islice(inf, 4))
			if not fastq_line_block:
				break
			if(len(fastq_line_block[1]) > max_length):
				max_length = len(fastq_line_block[1])
				max_line_length = fastq_line_block[1]
	return max_length

#change files to linux format
def change_to_linux_format(infile, outfile):
	with open(outfile, 'w') as outf:
		with open(infile, 'r') as inf:
			for line in inf:
				line.replace('\r\n', '\n')
				outf.write(line)

def copy_number_estimator_regression(rDNA_first_coord, rDNA_last_coord):
	#create temporary directory
	temp_dir = 'temp'+str(random.randint(10000, 100000))
	subprocess.Popen(['mkdir', temp_dir], stderr=subprocess.PIPE)	#creating the temporRY DIRECTORY
	subprocess.Popen(['cd', temp_dir], stderr=subprocess.PIPE)	#accessing the directory first before writing file in it to avoid error
	#change files to linux forma
	rDNA_fasta = temp_dir + '/rDNA_ref.fasta'
	change_to_linux_format(rDNA_ref, rDNA_fasta)
	whole_genome_fasta = temp_dir + '/whole_genome_ref.fasta'
	change_to_linux_format(whole_genome_ref, whole_genome_fasta)
	#find maximum read length
	max_read_length_fwd = max_line_length(fwd_reads_file)   #maximum read length in forward read file
	max_read_length_rvs = max_line_length(rvs_reads_file)   #maximum read length in reverse read file
	max_read_length = max(max_read_length_fwd, max_read_length_rvs) #greater of forward and reverse reads
	shoulder_size = max_read_length - 1
	rDNA_first_coord = int(rDNA_first_coord) - shoulder_size
	rDNA_last_coord = int(rDNA_last_coord) + shoulder_size
	#trim the rDNA fasta file
	rDNA_trm_fasta = temp_dir + '/trm_' + 'rDNA_ref.fasta'
	trim_fasta(rDNA_fasta, rDNA_trm_fasta, int('{0}'.format(rDNA_first_coord)), int('{0}'.format(rDNA_last_coord)))
	#index reference files
	try:
		subprocess.Popen(['bwa', 'index', rDNA_trm_fasta], stderr=subprocess.PIPE)
	except OSError:
		print('Error37: Cannot find bwa.')
		sys.exit(2)
	subprocess.Popen(['bwa', 'index', whole_genome_fasta], stderr=subprocess.PIPE).wait()
	#map the reads
	devnull = open(os.devnull, 'w')
	rDNA_sam = temp_dir + '/rDNA.sam'
	with open (rDNA_sam, 'w') as rDNA_samfile:
		subprocess.Popen(['bwa', 'mem', rDNA_trm_fasta, fwd_reads_file, rvs_reads_file], stdout=rDNA_samfile, stderr=devnull)	#to suppress screen output
	whole_genome_sam =  temp_dir + '/whole_genome.sam'
	with open (whole_genome_sam, 'w') as whole_genome_samfile:
		subprocess.Popen(['bwa', 'mem', whole_genome_fasta, fwd_reads_file, rvs_reads_file], stdout=whole_genome_samfile, stderr=devnull).wait()
	#sam to bam
	rDNA_bam = temp_dir + '/rDNA.bam'
	with open (rDNA_bam, 'w') as rDNA_bamfile:
		try:
			subprocess.Popen(['samtools', 'view', '-bS', rDNA_sam], stdout=rDNA_bamfile, stderr=subprocess.PIPE)
		except OSError:
			print('Error38: Cannot find samtools.')
			sys.exit(2)
	whole_genome_bam =  temp_dir + '/whole_genome.bam'
	with open (whole_genome_bam, 'w') as whole_genome_bamfile:
		subprocess.Popen(['samtools', 'view', '-bS', whole_genome_sam], stdout=whole_genome_bamfile, stderr=subprocess.PIPE).wait()
	#sort bam
	rDNA_sort_bam = temp_dir + '/rDNA_sort.bam'
	with open (rDNA_sort_bam, 'w') as rDNA_bamsort:
		subprocess.Popen(['samtools', 'sort', rDNA_bam], stdout=rDNA_bamsort, stderr=subprocess.PIPE)
	whole_genome_sort_bam = temp_dir + '/whole_genome_sort.bam'
	with open (whole_genome_sort_bam, 'w') as whole_genome_bamsort:
		subprocess.Popen(['samtools', 'sort', whole_genome_bam], stdout=whole_genome_bamsort, stderr=subprocess.PIPE).wait()
	#compute depth profile
	rDNA_depth = temp_dir + '/rDNA_depth.txt'
	with open (rDNA_depth, 'w') as rDNA_depthfile:
		subprocess.Popen(['samtools', 'depth', '-aa', rDNA_sort_bam], stdout=rDNA_depthfile, stderr=subprocess.PIPE)
	whole_genome_depth = temp_dir + '/whole_genome_depth.txt'
	with open (whole_genome_depth, 'w') as whole_genome_depthfile:
		subprocess.Popen(['samtools', 'depth', '-aa', whole_genome_sort_bam], stdout=whole_genome_depthfile, stderr=subprocess.PIPE).wait()
	#copy number estimate
	copy_number_estimate = examine(whole_genome_depth, int('{0}'.format(rDNA_start_locus)), int('{0}'.format(rDNA_end_locus)), int('{0}'.format(rDNA_buffer)), rDNA_depth, int('{0}'.format(shoulder_size)), int('{0}'.format(window_size)), int('{0}'.format(step_size)))
	print('copy_number_estimate = {0}'.format(copy_number_estimate))
	subprocess.Popen(['rm', '-rf', temp_dir], stderr=subprocess.PIPE)
	subprocess.Popen(['rm', '-rf', '__pycache__'], stderr=subprocess.PIPE)

#median based copy number estimator
def copy_number_estimator_median(rDNA_chromosome):
	#create temporary directory
	temp_dir = 'temp'+str(random.randint(10000, 100000))
	subprocess.Popen(['mkdir', temp_dir], stderr=subprocess.PIPE)   #creating the temporRY DIRECTORY
	subprocess.Popen(['cd', temp_dir], stderr=subprocess.PIPE)      #accessing the directory first before writing file in it to avoid error
	#change files to linux format
	whole_genome_fasta = temp_dir + '/whole_genome_ref.fasta'
	change_to_linux_format(whole_genome_ref, whole_genome_fasta)
	#index reference files
	try:
		 subprocess.Popen(['bwa', 'index', whole_genome_fasta], stderr=subprocess.PIPE).wait()
	except OSError:
		print('Error37: Cannot find bwa.')
		sys.exit(2)
	#map the reads
	devnull = open(os.devnull, 'w')
	whole_genome_sam =  temp_dir + '/whole_genome.sam'
	with open (whole_genome_sam, 'w') as whole_genome_samfile:
		subprocess.Popen(['bwa', 'mem', whole_genome_fasta, fwd_reads_file, rvs_reads_file], stdout=whole_genome_samfile, stderr=devnull).wait()
	#sam to bam
	whole_genome_bam =  temp_dir + '/whole_genome.bam'
	with open (whole_genome_bam, 'w') as whole_genome_bamfile:
		try:
			subprocess.Popen(['samtools', 'view', '-bS', whole_genome_sam], stdout=whole_genome_bamfile, stderr=subprocess.PIPE).wait()
		except OSError:
			print('Error38: Cannot find samtools.')
			sys.exit(2)
	whole_genome_sort_bam = temp_dir + '/whole_genome_sort.bam'
	with open (whole_genome_sort_bam, 'w') as whole_genome_bamsort:
		subprocess.Popen(['samtools', 'sort', whole_genome_bam], stdout=whole_genome_bamsort, stderr=subprocess.PIPE).wait()
	whole_genome_depth = temp_dir + '/whole_genome_depth.txt'
	with open (whole_genome_depth, 'w') as whole_genome_depthfile:
		subprocess.Popen(['samtools', 'depth', '-aa', whole_genome_sort_bam], stdout=whole_genome_depthfile, stderr=subprocess.PIPE).wait()
	#retain rDNA chromosome
	chrome_depth = temp_dir + '/chrome_depth.txt'
	rDNA_chromosome = '$1 == "' + rDNA_chromosome + '" { print }'
	with open (chrome_depth, 'w') as chrome_depthfile:
		subprocess.Popen(['awk', rDNA_chromosome, whole_genome_depth], stdout=chrome_depthfile, stderr=subprocess.PIPE).wait()
	#retain only rDNA region
	rDNA_depth = temp_dir + '/rDNA_depth.txt'
	rDNA_region = 'FNR >= ' + rDNA_start_locus + ' && FNR <= ' + rDNA_end_locus
	with open (rDNA_depth, 'w') as rDNA_depthfile:
		subprocess.Popen(['awk', rDNA_region, chrome_depth], stdout=rDNA_depthfile, stderr=subprocess.PIPE)
	#retain only non rDNA region
	nrDNA_depth = temp_dir + '/nrDNA_depth.txt'
	nrDNA_region = 'FNR < ' + rDNA_start_locus + ' || FNR > ' + rDNA_end_locus
	with open (nrDNA_depth, 'w') as nrDNA_depthfile:
		subprocess.Popen(['awk', nrDNA_region, chrome_depth], stdout=nrDNA_depthfile, stderr=subprocess.PIPE).wait()
	#sort rDNA depth file
	rDNA_depth_sorted = temp_dir + '/rDNA_depth_sorted.txt'
	with open (rDNA_depth_sorted, 'w') as rDNA_sortedfile:
		subprocess.Popen(['sort', '-nk', '3', rDNA_depth], stdout=rDNA_sortedfile, stderr=subprocess.PIPE)
	#sort nrDNA depth file
	nrDNA_depth_sorted = temp_dir + '/nrDNA_depth_sorted.txt'
	with open (nrDNA_depth_sorted, 'w') as nrDNA_sortedfile:
		subprocess.Popen(['sort', '-nk', '3', nrDNA_depth], stdout=nrDNA_sortedfile, stderr=subprocess.PIPE).wait()
	median = '{arr[NR]=$3} END {if (NR%2==1) print arr[(NR+1)/2]; else print (arr[NR/2]+arr[NR/2+1])/2}'
	rDNA_median = subprocess.Popen(['awk', median, rDNA_depth_sorted], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	rDNA_median = rDNA_median.stdout.readlines()
	rDNA_median = int(rDNA_median[0])
	nrDNA_median = subprocess.Popen(['awk', median, nrDNA_depth_sorted], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	nrDNA_median = nrDNA_median.stdout.readlines()
	nrDNA_median = int(nrDNA_median[0])
	copy_number_estimate = int(round(rDNA_median/nrDNA_median))
	print('copy_number_estimate = {0}'.format(copy_number_estimate))
	subprocess.Popen(['rm', '-rf', temp_dir], stderr=subprocess.PIPE)

#main function
def main():
	if reg:
		copy_number_estimator_regression(rDNA_first_coord, rDNA_last_coord)
	else:
		copy_number_estimator_median(rDNA_chromosome)
	
	

#function to verify input
def input_check(atribute, value, error_nr, digit):
	if value.startswith('-'):
		print('Error{0}: Value of {1} is not given.'.format(error_nr, atribute))
		sys.exit(2)
	if digit:
		if not value.isdigit():
			print('Error{0}: {1} must be a number.'.format(error_nr+1, atribute))
			sys.exit(2)
	return value

#class to take command line arguments
class estimate_copy_number():
	def __init__(self, args, args_len):
		global whole_genome_ref, rDNA_start_locus, rDNA_end_locus, rDNA_chromosome, rDNA_buffer, rDNA_ref, rDNA_first_coord, rDNA_last_coord, fwd_reads_file, rvs_reads_file, window_size, step_size, reg, med
		reg = False		#regression
		med = False		#median
		whole_ref = False	#whole genome reference
		rDNA_start = False	#start locus of rDNA in whole genome file
		rDNA_end = False	#end locus of rDNA in whole genome file
		rDNA_name = False	#name of rDNA chromosome
		bufer = False		#buffer between rDNA and non rDNA region
		ref = False		#reference file containing only rDNA
		first_coord = False	#first coordinate of rDNA unit
		last_coord = False	#last coordinate of rDNA unit
		fwd_reads = False	#forward reads file
		rvs_reads = False	#reverse reads file
		window = False		#size of window used in regression
		step = False		#size of step used in regression
		error = True		#There is error in input by default
		for i in range(args_len-1):		#loop of arguments
			if args[i].startswith('-'):	#this is an argument attribute
				if args[i] == '--regression' or args[i] == '-R':
					reg = True
				elif args[i] == '--median' or args[i] == '-M':
					med = True
				elif args[i] == '--whole_reference' or args[i] == '-F':
					whole_genome_ref = input_check('Whole genome reference file address', args[i+1], 8, False) 
					whole_ref = True
				elif args[i] == '--rDNA_start' or args[i] == '-U':
					rDNA_start_locus = input_check('rDNA start locus', args[i+1], 1, True)     #value will only be assigned after input check
					rDNA_start = True
				elif args[i] == '--rDNA_end' or args[i] == '-V':
					rDNA_end_locus = input_check('rDNA end locus', args[i+1], 3, True)     #value will only be assigned after input check
					rDNA_end = True
				elif args[i] == '--rDNA_name' or args[i] == '-n':
					rDNA_chromosome = input_check('rDNA chromosome name', args[i+1], 5, False)     #value will only be assigned after input check
					rDNA_name = True
				elif args[i] == '--buffer' or args[i] == '-b':
					rDNA_buffer = input_check('Buffer', args[i+1], 6, True)     #value will only be assigned after input check
					bufer = True
				elif args[i] == '--rDNA_reference' or args[i] == '-f':
					rDNA_ref = input_check('rDNA reference file address', args[i+1], 8, False)      #value will only be assigned after input check
					ref = True
				elif args[i] == '--unit_start' or args[i] == '-u':
					rDNA_first_coord = input_check('First coordinate of rDNA unit', args[i+1], 9, True)     #value will only be assigned after input check
					first_coord = True
				elif args[i] == '--unit_end' or args[i] == '-v':
					rDNA_last_coord = input_check('Last coordinate of rDNA unit', args[i+1], 11, True)     #value will only be assigned after input check
					last_coord = True
				elif args[i] == '--reads1' or args[i] == '-r':
					fwd_reads_file = input_check('First read file address', args[i+1], 13, False)      #value will only be assigned after input check
					fwd_reads = True
				elif args[i] == '--reads2' or args[i] == '-q':
					rvs_reads_file = input_check('Reverse read file address', args[i+1], 14, False)      #value will only be assigned after input check
					rvs_reads = True
				elif args[i] == '--window_size' or args[i] == '-w':
					window_size = input_check('Window size for regression', args[i+1], 15, True)     #value will only be assigned after input check
					window = True
				elif args[i] == '--step_size' or args[i] == '-s':
					step_size = input_check('Step size for regression', args[i+1], 17, True)     #value will only be assigned after input check
					step = True
				else:
					print('Error18: option {0} is not recognised.'.format(args[i]))
					sys.exit(2)
		#error list
		if reg and med:
			print('Error19: Only one of the options from regression or median should be given.')
		if not reg and not med:
			print('Error20: One of the options from regression or median must be given.')
		elif not whole_ref:
			print('Error21: Whole genome reference file is missing.')
		elif not rDNA_start:
			print('Error22: Start locus of rDNA is missing.')
		elif not rDNA_end:
			print('Error23: End locus of rDNA is missing.')
		elif med and not rDNA_name:
			print('Error24: rDNA chromosome name is missing.')
		elif reg and not bufer:
			print('Error25: Buffer between rDNA and non rDNA region is missing.')
		elif reg and not ref:
			print('Error26: rDNA reference file is missing.')
		elif reg and not first_coord:
			print('Error27: First coordinate of rDNA unit is missing.')
		elif reg and not last_coord:
			print('Error28: Last coordinate of rDNA unit is missing.')
		elif not fwd_reads:
			print('Error29: Read1 file is missing.')
		elif not rvs_reads:
			print('Error30: Read2 file is missing.')
		elif reg and not window:
			print('Error31: Window size for regression is missing.')
		elif reg and not step:
			print('Error32: Step size for regression is missing.')
		elif not path.exists(whole_genome_ref):
			print('Error33: Cannot find whole genome reference file.')
		elif reg and not path.exists(rDNA_ref):
			print('Error34: Cannot find rDNA reference file.')
		elif not path.exists(fwd_reads_file):
			print('Error35: Cannot find reads1 file.')
		elif not path.exists(rvs_reads_file):
			print('Error36: Cannot find reads2 file.')
		else:
			error = False
		if error:
			sys.exit(2)
	def execute(self):
		main()


if __name__ == "__main__":
	full_cmd_arguments = sys.argv   # Get full command-line arguments
	argument_list = full_cmd_arguments[1:]  # Keep all but the first (because first is the function name)
	short_options = "R:M:F:U:V:n:b:f:u:v:r:q:w:s:"
	long_options = ["regression=", "median=", "whole_reference=", "rDNA_start=", "rDNA_end=", "rDNA_name=", "buffer=", "rDNA_reference=", "unit_start=", "unit_end=", "reads1=", "reads2=", "window_size=", "step_size="]
	try:
		arguments, values = getopt.getopt(argument_list, short_options, long_options)
	except getopt.error as err:
		# Output error, and return with an error code
		print (str(err))
		sys.exit(2)
	args_len = len(sys.argv)
	#args = []
	#for i in range(1, args_len):
	#	args.append(sys.argv[i])
	estimate_copy_number(sys.argv[1:], args_len).execute()

