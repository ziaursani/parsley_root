import sys, os, subprocess, random
from os import path
from itertools import islice
from fasta_trimmer import trim_fasta
from copyno import examine

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

#maximum read length
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

#estimate copy number
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
		exit_program('Error9: Cannot find bwa.')
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
			exit_program('Error9: Cannot find samtools.')
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

#main function
def main():
	copy_number_estimator_regression(rDNA_first_coord, rDNA_last_coord)
		
#class to take command line arguments
class estimate_copy_number_regression():
	def __init__(self, args, args_len):
		global whole_genome_ref, rDNA_start_locus, rDNA_end_locus, rDNA_buffer, rDNA_ref, rDNA_first_coord, rDNA_last_coord, fwd_reads_file, rvs_reads_file, window_size, step_size
		whole_ref = False	#whole genome reference
		rDNA_start = False	#start locus of rDNA in whole genome file
		rDNA_end = False	#end locus of rDNA in whole genome file
		bufer = False		#buffer between rDNA and non rDNA region
		ref = False		#reference file containing only rDNA
		first_coord = False	#first coordinate of rDNA unit
		last_coord = False	#last coordinate of rDNA unit
		fwd_reads = False	#forward reads file
		rvs_reads = False	#reverse reads file
		window = False		#size of window used in regression
		step = False		#size of step used in regression
		for i in range(args_len-1):		#loop of arguments
			if args[i].startswith('-'):	#this is an argument attribute
				if args[i] == '--whole_reference' or args[i] == '-F':
					whole_genome_ref, whole_ref = input_check(args, i, 'str') 
				elif args[i] == '--rDNA_start' or args[i] == '-U':
					rDNA_start_locus, rDNA_start = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--rDNA_end' or args[i] == '-V':
					rDNA_end_locus, rDNA_end = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--buffer' or args[i] == '-b':
					rDNA_buffer, bufer = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--rDNA_reference' or args[i] == '-f':
					rDNA_ref, ref = input_check(args, i, 'str')      #value will only be assigned after input check
				elif args[i] == '--unit_start' or args[i] == '-u':
					rDNA_first_coord, first_coord = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--unit_end' or args[i] == '-v':
					rDNA_last_coord, last_coord = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--reads1' or args[i] == '-r':
					fwd_reads_file, fwd_reads = input_check(args, i, 'str')      #value will only be assigned after input check
				elif args[i] == '--reads2' or args[i] == '-q':
					rvs_reads_file, rvs_reads = input_check(args, i, 'str')      #value will only be assigned after input check
				elif args[i] == '--window_size' or args[i] == '-w':
					window_size, window = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--step_size' or args[i] == '-s':
					step_size, step = input_check(args, i, 'int')     #value will only be assigned after input check
				else:
					exit_program('Error0: Argument not recognised '+args[i])

		if not whole_ref:
			exit_program('Error6: Missing argument: --whole_reference or -F.')
		elif not rDNA_start:
			exit_program('Error6: Missing argument: --rDNA_start or -U.')
		elif not rDNA_end:
			exit_program('Error6: Missing argument: --rDNA_end or -V.')
		elif not bufer:
			exit_program('Error6: Missing argument: --buffer or -b.')
		elif not ref:
			exit_program('Error6: Missing argument: --rDNA_reference or -f.')
		elif not first_coord:
			exit_program('Error6: Missing argument: --unit_start or -u.')
		elif not last_coord:
			exit_program('Error6: Missing argument: --unit_end or -v.')
		elif not fwd_reads:
			exit_program('Error6: Missing argument: --reads1 or -r')
		elif not rvs_reads:
			exit_program('Error6: Missing argument: --reads2 or -q')
		elif not window:
			exit_program('Error6: Missing argument: --window_size or -w.')
		elif not step:
			exit_program('Error6: Missing argument: --step_size or -s.')
			
	def execute(self):
		main()
