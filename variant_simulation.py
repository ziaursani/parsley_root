import sys, os, subprocess, configparser, random
from os import path
from fasta_trimmer import trim_fasta
from genome_generator import generate_genome
from variant_simulator import generate_variants

#change files to linux format
def change_to_linux_format(infile, outfile):
        with open(outfile, 'w') as outf:
                with open(infile, 'r') as inf:
                        for line in inf:
                                line.replace('\r\n', '\n')
                                outf.write(line)

#function to exit with error message
def exit_program(error_message):
        print(error_message)
        sys.exit(2)

#function to verify input
def input_check(arg_name, arg_value, Type):
        if Type == 'int':
                if not arg_value.isdigit():
                        exit_program('Error2: Not a number: '+arg_name)
        elif Type == 'float':
                if not arg_value.replace('.','',1).isdigit():
                        exit_program('Error3: Not numeric: '+arg_name)
                if float(arg_value) < 0.0 or float(arg_value) > 1.0:
                        exit_program('Error4: Out of range: '+arg_name)
        elif Type == 'str':
                if not path.exists(arg_value):
                        exit_program('Error5: Does not exist: '+arg_name)
        return arg_value

def reads_generator(info):
	rDNA_structure = info['rDNA structure']
	nr_of_repeats = input_check('Number of rDNA unit repeats', rDNA_structure['Number of rDNA unit repeats'], 'int')
	input_reference_file = info['Input reference file']
	ref_file = input_check('Name and address of input reference fasta file', input_reference_file['Name and address of input reference fasta file'], 'str')
	first_coordinate = input_check('First Coordinate of rDNA unit', input_reference_file['First Coordinate of rDNA unit'], 'int')
	last_coordinate = input_check('Last Coordinate of rDNA unit', input_reference_file['Last Coordinate of rDNA unit'], 'int')
	output_svcf_genome_files = info['Output simulated vcf and genome files']
	simulated_vcf_file = output_svcf_genome_files['Name and address of output simulated vcf file']
	simulated_rDNA_file = output_svcf_genome_files['Name and address of output simulated rDNA fasta file']
	length_of_line = input_check('Length of line of rDNA file', output_svcf_genome_files['Length of line of rDNA file'], 'int')
	Output_read_files = info['Output read files']
	read_depth = input_check('Read Depth', Output_read_files['Read Depth'], 'int')
	snp_error_rate = input_check('SNP Error Rate', Output_read_files['SNP Error Rate'], 'float')
	read_file_prefix = Output_read_files['Read Filename Prefix']
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
	#generating the simulated vcf
	generate_variants(info, temp_dir).execute()
	#generating the simulated rDNA
	generate_genome(temp_dir, '{0}'.format(simulated_vcf_file), '{0}'.format(simulated_rDNA_file), int('{0}'.format(nr_of_repeats)), int('{0}'.format(length_of_line)))
	#activating the pIRS
	command = "export LD_LIBRARY_PATH=/usr/local/lib64:${LD_LIBRARY_PATH}"
	subprocess.Popen([command], shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	#generating reads
	devnull = open(os.devnull, 'w')	#to suppress screen output
	try:
		subprocess.Popen(['pirs', 'simulate', '{0}'.format(simulated_rDNA_file), '-x', '{0}'.format(read_depth), '-e', '{0}'.format(snp_error_rate), '-o', '{0}'.format(read_file_prefix)], stdout=devnull, stderr=devnull)
	except OSError:
		exit_program('Error9: Cannot find: pirs.')
	subprocess.Popen(['rm', '-rf', temp_dir], stderr=subprocess.STDOUT)	
	
#main function
def simulate_variants(config):
	config = input_check('configuration file', config, 'str')
	info = configparser.ConfigParser()
	info.read(config)
	reads_generator(info)
