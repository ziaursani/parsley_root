import sys, os, subprocess, random
from os import path
from itertools import islice
from fasta_trimmer import trim_fasta

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
		exit_program('Error9: Cannot find bwa.')
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
			exit_program('Error9: Cannot find samtools.')
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

def main():
	copy_number_estimator_median(rDNA_chromosome)

#class to take command line arguments
class estimate_copy_number_median():
	def __init__(self, args, args_len):
		global whole_genome_ref, rDNA_start_locus, rDNA_end_locus, rDNA_chromosome, fwd_reads_file, rvs_reads_file
		whole_ref = False	#whole genome reference
		rDNA_start = False	#start locus of rDNA in whole genome file
		rDNA_end = False	#end locus of rDNA in whole genome file
		rDNA_name = False	#name of rDNA chromosome
		ref = False		#reference file containing only rDNA
		fwd_reads = False	#forward reads file
		rvs_reads = False	#reverse reads file
		for i in range(args_len-1):		#loop of arguments
			if args[i].startswith('-'):	#this is an argument attribute
				if args[i] == '--whole_reference' or args[i] == '-F':
					whole_genome_ref, whole_ref = input_check(args, i, 'str')
				elif args[i] == '--rDNA_start' or args[i] == '-U':
					 rDNA_start_locus, rDNA_start = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--rDNA_end' or args[i] == '-V':
					rDNA_end_locus, rDNA_end = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--rDNA_name' or args[i] == '-n':
					rDNA_chromosome, rDNA_name = input_check(args, i, 'str*')
				elif args[i] == '--reads1' or args[i] == '-r':
					fwd_reads_file, fwd_reads = input_check(args, i, 'str')      #value will only be assigned after input check
				elif args[i] == '--reads2' or args[i] == '-q':
					rvs_reads_file, rvs_reads = input_check(args, i, 'str')      #value will only be assigned after input check
				else:
					exit_program('Error0: Argument not recognised '+args[i])
		#error list
		if not whole_ref:
			exit_program('Error6: Missing argument: --whole_reference or -F.')
		elif not rDNA_start:
			exit_program('Error6: Missing argument: --rDNA_start or -U.')
		elif not rDNA_end:
			exit_program('Error6: Missing argument: --rDNA_end or -V.')
		elif not rDNA_name:
			exit_program('Error6: Missing argument: --rDNA_name or -n.')
		elif not fwd_reads:
			exit_program('Error6: Missing argument: --reads1 or -r.')
		elif not rvs_reads:
			exit_program('Error6: Missing argument: --reads2 or -q.')
	def execute(self):
		main()
