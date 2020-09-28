# parsley_root
This is pre-release.
##
Parsley Root: Pipeline for Analysis of Ribosomal Locus Evolution in Yeast (Reusable for Other Organisms Though)
Version: V0
Release Date: September 2020
Dependencies: Python, BWA, Samtools including Bcftools, Freebayes, Vcflib, pIRS, Samclip, FAT-Cigar, Tabix
Instructions for use
Download from top right corner
Unzip the folder
Create link
chmod +x parsley.py
ln -s parsley.py parsley
Type ./parsley --help
or follow instructions below:
These commands are mutually exclusive
        discover_variants       Discovers variants in strain against reference sequence.
        simulate_variants       Simulates variant strains from reference sequence.
        compare_variants        Compares variants in simulated strain against predicted variants of the strain.
        estimate_copy_number    Estimates number of copies of Ribosomal DNA unit
Usage:  ./parsley discover_variants [arguments]
About:  This discovers variants in the Illumina reads from variant strain against a reference sequence of Ribosomal DNA. It gives output in vcf format.
Arguments:
        Mandatory:
                -f      --rDNA_reference        Reference fasta file of rDNA
                -u      --unit_start            Start locus of rDNA unit in the reference file
                -v      --unit_end              End locus of rDNA unit in the reference file
                -r      --reads1                Forward reads fastq file
                -q      --reads2                Reverse reads fastq file
                -o      --out                   Output VCF file
        Optional:
                -s      --sub_sample            Takes (0 < x < 1) portion of reads at random
Usage:  ./parsley simulate_variants [arguments]
About:  This outputs genome, reads, and variant list of a simulated variant.
Arguments:
 Mandatory:
                -       -                       Configuration file present in the parsley repository. Change parameters as they suite the strains under your study.
Usage:   ./parsley compare_variants [arguments]
About:   This compares variants in the simulated strain against predicted variants, by reporting false positives and true negatives.
Arguments:
 Mandatory:
                -I      --simulated_vcf         variant list of simulated strain obtained from "simulate_variants" command
                -i      --predictive_vcf        variant list of simulated strain obtained from "discover_variants" command
                -o      --comparative_vcf       comparison of variants of other two input variant lists (output)
                -n      --copy_number           number of copies of rDNA unit in the simulated variant
Usage:  ./parsley estimate_copy_number [options] [arguments]
About:  This outputs number of copies of rDNA unit.
Options:        These options are mutually exclusive
        -M      --median        Statistical median of depth over the rDNA and non-rDNA regions of the chromosome
        -R      --regression    Estimating number of copies based on regression analysis
Usage:  ./parsley estimate_copy_number --median [arguments]
About:  This outputs number of copies of rDNA unit based on median method
Arguments:
        Mandatory:
                -F      --whole_reference       Reference fasta file of whole genome
                -U      --rDNA_start            Start locus of rDNA in the chromosome
                -V      --rDNA_end              End locus of rDNA in the chromosome
                -n      --rDNA_name             Name of chromosome containing rDNA
                -r      --reads1                Farward reads fastq file
                -q      --reads2                Reverse reads fastq file
Usage:   ./parsley estimate_copy_number --regression [arguments]
About:  This outputs number of copies of rDNA unit based on regression method
Arguments:
        Mandatory:
        -F      --whole_reference       Reference fasta file of whole genome
        -U      --rDNA_start            Start locus of rDNA in the chromosome
        -V      --rDNA_end              End locus of rDNA in the chromosome
        -b      --buffer                Number of nucleotides between rDNA and non-rDNA region
        -f      --rDNA_referene         Reference fasta file for rDNA
        -u      --unit_start            Starting nucleotide of rDNA unit
        -v      --unit_end              End nucleotide of rDNA unit
        -r      --reads1                Farward reads fastq file
        -q      --reads2                Reverse reads fastq file
        -w      --window_size           Window size for sliding function
        -s      --step_size             Step size for sliding function
