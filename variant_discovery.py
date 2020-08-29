import getopt, sys
import subprocess, random
from os import path
from itertools import islice
from fasta_trimmer import trim_fasta
from variant_unifier import unify_variants

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

def variant_discoverer(first_coordinate, last_coordinate, max_read_len, rDNA_unit_size):
	#create temporary directory
	temp_dir = 'temp'+str(random.randint(10000, 100000))
	subprocess.Popen(['mkdir', temp_dir], stderr=subprocess.PIPE)   #creating the temporRY DIRECTORY
	subprocess.Popen(['cd', temp_dir], stderr=subprocess.PIPE)      #accessing the directory first before writing file in it to avoid error
	#change file to linux format
	linux_ref = temp_dir + '/linux.fasta'
	change_to_linux_format(ref_file, linux_ref)
	#trim reference file
	trim_ref = temp_dir + '/trim.fasta'
	trim_fasta(linux_ref, trim_ref, int('{0}'.format(first_coordinate)), int('{0}'.format(last_coordinate)))
	#Index fasta with bwa
	try:
		subprocess.Popen(['bwa', 'index', trim_ref], stderr=subprocess.PIPE).wait()
	except OSError:
		exit_program('Error9: Cannot find bwa.')
	#Index fasta with samtool
	try:
		subprocess.Popen(['samtools', 'faidx', trim_ref], stderr=subprocess.PIPE)
	except OSError:
		exit_program('Error9: Cannot find samtools.')
	#Map the reads with bwa
	process_align_reads = subprocess.Popen(['bwa', 'mem', trim_ref,  fwd_reads_file,  rvs_reads_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#filter unmapped reads with samtools
	process_filter_unmapped_reads = subprocess.Popen(['samtools', 'view', '-h', '-F', '4'], stdin=process_align_reads.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#filter reads with soft clippings with software available at https://github.com/tseemann/samclip
	try:
		process_filter_softclipped_reads = subprocess.Popen(['samclip', '--ref', trim_ref, '--max', '0'], stdin=process_filter_unmapped_reads.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		exit_program('Error9: Cannot find samclip.')
	#Convert file to BAM format for the next tool
	bamfile = temp_dir + '/bamfile.bam'
	with open (bamfile, 'w') as bam:
		subprocess.Popen(['samtools', 'view', '-bS'], stdin=process_filter_softclipped_reads.stdout, stdout=bam, stderr=subprocess.PIPE).wait()
	#Use Prithika's new Python tool to convert CIGAR strings to CIGAR2 strings software available from Prithika
	fat_cigar_bam = temp_dir + '/fat_cigar.bam'
	subprocess.Popen([python, 'fat-cigar.py', 'linear', bamfile, fat_cigar_bam], stderr=subprocess.PIPE).wait()
	#Convert the BAM files back to SAM format for CIGAR2 filtering
	fat_cigar_sam = temp_dir + '/fat_cigar.sam'
	with open (fat_cigar_sam, 'w') as samfile:
		subprocess.Popen(['samtools', 'view', '-h', fat_cigar_bam], stdout=samfile, stderr=subprocess.PIPE).wait()
	#Copy header from the sam file
	copy_header = temp_dir + '/ch.sam'
	with open (copy_header, 'w') as chfile:
		subprocess.Popen(['samtools', 'view', '-H', fat_cigar_sam], stdout=chfile, stderr=subprocess.PIPE)
	#CIGAR2 filtering - only keep reads if they they have at least 20 nucleotide matches at each end of the read
	filtered_sam = temp_dir + '/filtered.sam'
	with open (filtered_sam, 'w') as filteredfile:
		subprocess.Popen(['awk', '($6 ~ /^[2-9][0-9]=.*[2-9][0-9]=$/) || ($6 ~ /^[2-9][0-9]=.*1[0-9]{2}=$/) || ($6 ~ /^1[0-9]{2}=.*[2-9][0-9]=$/) || ($6 ~ /^[1-9][0-9]{1,2}=$/)', fat_cigar_sam], stdout=filteredfile, stderr=subprocess.PIPE).wait()
	#Add the header to the CIGAR2-filtered files as the filtering step lost it
	process_add_header = subprocess.Popen(['cat', copy_header, filtered_sam], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#Convert the CIGAR2-filtered SAM files to BAM format
	back_to_bam = temp_dir + '/back.bam'
	with open (back_to_bam, 'w') as bam:
		if float(sub_sample) < 1.0:
			subprocess.Popen(['samtools', 'view', '-bS', '-s', sub_sample], stdin=process_add_header.stdout, stdout=bam, stderr=subprocess.PIPE).wait()
		else:
			subprocess.Popen(['samtools', 'view', '-bS'], stdin=process_add_header.stdout, stdout=bam, stderr=subprocess.PIPE).wait()
	#Stack the filtered reads under the reference - up to maximum depth using freebayes
	vcf_file = temp_dir + '/vcfile.vcf'
	with open (vcf_file, 'w') as vcfile:
		try:
			subprocess.Popen(['freebayes', '-f', trim_ref, '-F', '0.001', '--pooled-continuous', back_to_bam], stdout=vcfile, stderr=subprocess.PIPE).wait()
		except OSError:
			exit_program('Error9: Cannot find freebayes.')
	#Break multi-allelic sites
	try:
		process_breakmulti = subprocess.Popen(['vcfbreakmulti', vcf_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		exit_program('Error9: Cannot find vcflib.')
	#Decompose biallelic sites
	process_allelicprimitives = subprocess.Popen(['vcfallelicprimitives', '-kg'], stdin=process_breakmulti.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#Fix the field of alternate allelic depth
	process_fix_depth = subprocess.Popen(['./vcf_ad_fix.sh'], stdin=process_allelicprimitives.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#normalise the vcf
	try:
		process_normalise = subprocess.Popen(['bcftools', 'norm', '-f', trim_ref, '-m-'], stdin=process_fix_depth.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	except OSError:
		exit_program('Error9: Cannot find bcftools.')
	#relocate variants
	relocated_vcf = temp_dir + '/relocated.vcf'
	subprocess.Popen([python, 'variant_relocator.py', relocated_vcf, str(max_read_len), str(rDNA_unit_size+max_read_len-1), '0'], stdin=process_normalise.stdout, stderr=subprocess.PIPE).wait()
	#sort variants
	sorted_vcf =  temp_dir + '/sorted.vcf'
	with open (sorted_vcf, 'w') as sortedfile:                      #file object to write on the file
		subprocess.Popen(['bcftools', 'sort', relocated_vcf], stdout=sortedfile, stderr=subprocess.PIPE).wait()
	#unify decomposed variants
	unified_dec_vcf = temp_dir + '/unified_dec.vcf'
	unify_variants(sorted_vcf, unified_dec_vcf, 1)
	#unify relocated variants
	unified_rel_vcf = temp_dir + '/unified_rel.vcf'
	unify_variants(unified_dec_vcf, unified_rel_vcf, 0)
	#filter variants
	filtered_vcf = temp_dir + '/filtered.vcf'
	with open (filtered_vcf, 'w') as filteredfile:                #file object to write on the file
		subprocess.Popen(['vcffilter', '-f', 'AO > 1 & AO / DP > 0.005', unified_rel_vcf], stdout=filteredfile, stderr=subprocess.PIPE).wait() #filtering the variants
	#retrieve frquency
	subprocess.Popen([python, 'frequency_retriever.py', temp_dir, str(max_read_len-1), str(rDNA_unit_size), out_vcf]).wait()
	subprocess.Popen(['rm', '-rf', temp_dir], stderr=subprocess.PIPE)
	
#main function
def main(first_coordinate, last_coordinate):
	max_read_length_fwd = max_line_length(fwd_reads_file)   #maximum read length in forward read file
	max_read_length_rvs = max_line_length(rvs_reads_file)   #maximum read length in reverse read file
	max_read_length = max(max_read_length_fwd, max_read_length_rvs) #greater of forward and reverse reads
	rDNA_unit_size = int(last_coordinate) - int(first_coordinate) + 1  #size of rDNA unit
	last_coordinate = int(last_coordinate) + max_read_length - 1	#last coordinate of trimmed fasta file
	first_coordinate = int(first_coordinate) - max_read_length + 1	#first coordinate of trimmed fasta file
	variant_discoverer(first_coordinate, last_coordinate, max_read_length,  rDNA_unit_size)	#fuction to discover the variants

#class to take command line arguments
class discover_variants():
	def __init__(self, args, args_len):
		global ref_file, first_coordinate, last_coordinate, fwd_reads_file, rvs_reads_file, sub_sample, out_vcf, python
		ref = False
		first = False
		last = False
		fwd = False
		rvs = False
		out = False
		sub = False
		sub_sample = '1.0'			#default value
		for i in range(args_len-1):
			if args[i].startswith('-'):     #this is an argument attribute
				if args[i] == '--rDNA_reference' or args[i] == '-f':
					ref_file, ref = input_check(args, i, 'str')			#reference file
				elif args[i] == '--unit_start' or args[i] == '-u':
					first_coordinate, first = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--unit_end' or args[i] == '-v':
					last_coordinate, last = input_check(args, i, 'int')     #value will only be assigned after input check
				elif args[i] == '--reads1' or args[i] == '-r':
					fwd_reads_file, fwd = input_check(args, i, 'str')                      #read1 file
				elif args[i] == '--reads2' or args[i] == '-q':
					rvs_reads_file, rvs = input_check(args, i, 'str')                      #read2 file
				elif args[i] == '--sub_sample' or args[i] == '-s':
					sub_sample, sub = input_check(args, i, 'float')
				elif args[i] == '--out' or args[i] == '-o':
					out_vcf, out = input_check(args, i, 'str')                      #output vcf file
				else:
					exit_program('Error0: Argument not recognised '+args[i])
		#error and warning list
		if not ref:
			exit_program('Error6: Missing argument: --rDNA_reference or -f.')
		if not first:
			exit_program('Error6: Missing argument: --unit_start or -u.')
		if not last:
			exit_program('Error6: Missing argument: --unit_end or -v.')
		if not fwd:
			exit_program('Error6: Missing argument: --reads1 or -r')
		if not rvs:
			exit_program('Error6: Missing argument: --reads2 or -q')
		if not out:
			exit_program('Error6: Missing argument: out or -o')
		python = 'python'+str(sys.version_info.major)   #python version
	def execute(self):
		main(first_coordinate, last_coordinate)		#adding two arguments because these are not acting as global variables for unknown reasons
