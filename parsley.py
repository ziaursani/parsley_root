#!/usr/bin/env python
import sys, getopt, subprocess, os
from os import path
from copy_number_estimation_median import estimate_copy_number_median
from copy_number_estimation_regression import estimate_copy_number_regression
from variant_discovery import discover_variants
from variant_simulation import simulate_variants
from variant_comparator import compare_variants

#function to exit with error message
def exit_program(error_message):
	print(error_message)
	print('Type:	./parsley --help for further instructions.')
	subprocess.Popen(['rm', '-rf', '__pycache__'], stderr=subprocess.PIPE)
	sys.exit(2)

#function to verify input
def input_check(attribute, value, Type):
	if value.startswith('-'):
		exit_program('Error1: Value missing: '+attribute)
	elif Type == 'int':
		if not value.isdigit():
			exit_program('Error2: Not a number: '+attribute)
	elif Type == 'float':
		if not value.replace('.','',1).isdigit():
			exit_program('Error3: Not numeric: '+attribute)
		if float(value) < 0.0 or float(value) > 1.0:
			exit_program('Error4: Out of range: '+attribute)
	elif Type == 'str':
		if not path.exists(value):
			exit_program('Error5: Does not exist: '+attribute)
	return value, True

#function to generate help information
def help_info(args):
	if len(args) == 2:
		print('Parsley Root: Pipeline for Analysis of Ribosomal Locus Evolution in Yeast (Reusable for Other Organisms Though)')
		print('Version: V0')
		print('Release Date: September 2020')
		print('Dependencies: Python, BWA, Samtools including Bcftools, Freebayes, Vcflib, pIRS, Samclip, FAT-Cigar, Tabix')
		print('Instructions for use')
		print('Download from top right corner')
		print('Unzip the folder')
		print('Start using straight away')
		print('Usage: ./parsley <command> [Arguments]')
		print('Commands:			These commands are mutually exclusive')
		print('	discover_variants	Discovers variants in strain against reference sequence.')
		print('	simulate_variants	Simulates variant strains from reference sequence.')
		print('	compare_variants	Compares variants in simulated strain against predicted variants of the strain.')
		print('	estimate_copy_number	Estimates number of copies of Ribosomal DNA unit')
		print('Type:	./parsley --help <Command> for further instructions.')
	elif args[2] == 'discover_variants':
		print('About:	This discovers variants in the Illumina reads from variant strain against a reference sequence of Ribosomal DNA. It gives output in vcf format.')
		print('Usage:	./parsley discover_variants [arguments]')
		print('Arguments:')
		print('	Mandatory:')
		print('		-f	--rDNA_reference	Reference fasta file of rDNA')
		print('		-u	--unit_start		Start locus of rDNA unit in the reference file')
		print('		-v	--unit_end		End locus of rDNA unit in the reference file')
		print('		-r	--reads1		Forward reads fastq file')
		print('		-q	--reads2		Reverse reads fastq file')
		print('		-o	--out			Output VCF file')
		print('	Optional:')
		print('		-s	--sub_sample		Takes (0 < x < 1) portion of reads at random')
	elif args[2] == 'simulate_variants':
		print('About:	This outputs genome, reads, and variant list of a simulated variant.')
		print('Usage:	./parsley simulate_variants [arguments]')
		print('Arguments:')
		print(' Mandatory:')
		print('		-	-			Configuration file present in the parsley repository. Change parameters as they suite the strains under your study.')
	elif args[2] == 'compare_variants':
		print('About:   This compares variants in the simulated strain against predicted variants, by reporting false positives and true negativss.')
		print('Usage:   ./parsley compare_variants [arguments]')
		print('Arguments:')
		print(' Mandatory:')
		print('		-I      --simulated_vcf		variant list of simulated strain obtained from "simulate_variants" command')
		print('		-i	--predictive_vcf	variant list of simulated strain obtained from "discover_variants" command')
		print('		-o      --comparative_vcf	comparison of variants of other two input variant lists (output)')
		print('		-n	--copy_number		number of copies of rDNA unit in the simulated variant')
	elif args[2] == 'estimate_copy_number':
		if len(args) == 3:
			print('About:	This outputs number of copies of rDNA unit.')
			print('Usage:	./parsley estimate_copy_number [options] [arguments]')
			print('Options:	These options are mutually exclusive')
			print('	-M	--median	Statistical median of depth over the rDNA and non-rDNA regions of the chromosome')
			print('	-R	--regression	Estimating number of copies based on regression analysis')
			print('Type	./parsley --help estimate_copy_number [option]	for further instructions.')
		else:
			if args[3] == '--median' or args[3] == '-M':
				print('About:	This outputs number of copies of rDNA unit based on median method')
				print('Usage:	./parsley estimate_copy_number --median [arguments]')
				print('Arguments:')
				print('	Mandatory:')
				print('		-F	--whole_reference	Reference fasta file of whole genome')
				print('		-U	--rDNA_start		Start locus of rDNA in the chromosome')
				print('         	-V      --rDNA_end		End locus of rDNA in the chromosome')
				print('		-n	--rDNA_name		Name of chromosome containing rDNA')
				print('		-r	--reads1		Farward reads fastq file')
				print('         	-q      --reads2                Reverse reads fastq file')
			elif args[3] == '--regression' or args[3] == '-R':
				print('About:	This outputs number of copies of rDNA unit based on regression method')
				print('Usage:   ./parsley estimate_copy_number --regression [arguments]')
				print('Arguments:')
				print('	Mandatory:')
				print('	-F      --whole_reference       Reference fasta file of whole genome')
				print('	-U      --rDNA_start            Start locus of rDNA in the chromosome')
				print('	-V      --rDNA_end              End locus of rDNA in the chromosome')
				print('	-b	--buffer		Number of nucleotides between rDNA and non-rDNA region')
				print('	-f	--rDNA_referene		Reference fasta file for rDNA')
				print('	-u	--unit_start		Starting nucleotide of rDNA unit')
				print('	-v	--unit_end		End nucleotide of rDNA unit')
				print('	-r	--reads1		Farward reads fastq file')
				print('	-q      --reads2                Reverse reads fastq file')
				print('	-w	--window_size		Window size for sliding function')
				print('	-s	--step_size		Step size for sliding function')
			else:
				exit_program('Error0: Option not recognised: '+sys.argv[3])
	else:
		exit_program('Error0: Command not recognised: '+sys.argv[2])

#remove junk files
def remove_junk():
        directory = os.getcwd()
        files_in_directory = os.listdir(directory)
        for junk in files_in_directory:
                if junk.endswith(".pyc"):
                        os.remove(junk)

#main function
def main():
	if sys.argv[1] == 'estimate_copy_number':
		if sys.argv[2] == '--median' or sys.argv[2] == '-M':
			estimate_copy_number_median(sys.argv[3:], len(sys.argv)-2).execute()
		else:
			estimate_copy_number_regression(sys.argv[3:], len(sys.argv)-2).execute()
	elif sys.argv[1] == 'discover_variants':
		discover_variants(sys.argv[2:], len(sys.argv)-1).execute()
	elif sys.argv[1] == 'simulate_variants':
		simulate_variants(sys.argv[2])
	elif sys.argv[1] == 'compare_variants':
		compare_variants(sys.argv[2:], len(sys.argv)-1).execute()
	else:
		help_info(sys.argv)
	subprocess.Popen(['rm', '-rf', '__pycache__'], stderr=subprocess.PIPE)
	remove_junk()

if __name__ == "__main__":
	if len(sys.argv) < 2:
		exit_program('Error7: No command given.')
	full_cmd_arguments = sys.argv   # Get full command-line arguments
	argument_list = full_cmd_arguments[1:]  # Keep all but the first (because first is the function name)
	if (sys.argv[1] == 'discover_variants'):
		short_options = "f:u:v:r:q:s:o:"
		long_options = ["reference=", "unit_start=", "unit_end=", "reads1=", "reads2=", "sub_sample=", "out="]
	elif (sys.argv[1] == 'estimate_copy_number'):
		if len(sys.argv) < 3:
			exit_program('Error9 Command without any option.')
		if (sys.argv[2] == '--regression' or sys.argv[2] == '-R'):
			short_options = "R:F:U:V:b:f:u:v:r:q:w:s:"
			long_options = ["regression=", "whole_reference=", "rDNA_start=", "rDNA_end=", "buffer=", "rDNA_reference=", "unit_start=", "unit_end=", "reads1=", "reads2=", "window_size=", "step_size="]
		elif (sys.argv[2] == '--median' or sys.argv[2] == '-M'):
			short_options = "M:F:U:V:n:r:q:"
			long_options = ["median=", "whole_reference=", "rDNA_start=", "rDNA_end=", "rDNA_name=", "reads1=", "reads2="]
		else:
			exit_program('Error0: Option not recognised: '+sys.argv[2])
	elif (sys.argv[1] == 'compare_variants'):
		short_options = "I:i:o:n"
		long_options = ["simulated_vcf=", "predictive_vcf=", "comparative_vcf=", "copy_number"]
	elif (sys.argv[1] == 'simulate_variants'):
		short_options = ""
		long_options = []
	elif (sys.argv[1] == '--help' or sys.argv[1] == '-h'):
		short_options = "h:M:R"
		long_options = ["help=", "median=", "regression="]
	else:
		exit_program('Error0: Command not recognised: '+sys.argv[1])
	if len(sys.argv) >= 3 and sys.argv[1] != '--help' and sys.argv[1] != '-h':
		try:
			arguments, values = getopt.getopt(argument_list, short_options, long_options)
		except getopt.error as err:
			# Output error, and return with an error code
			exit_program(str(err))
	elif (not sys.argv[1] == '--help' and not sys.argv[1] == '-h'):
		exit_program('Error8: Command without any argument.')
	main()
