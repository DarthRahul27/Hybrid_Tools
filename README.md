# Hybrid_Tools
This is a set of perl and R scripts which allows us to analyse Hybrid Species 

## VCF-Public-Private-Alleles-Script.pl
This script uses VCF files to get single nucleotide polymorphism data, assigning those snps to designated lineage and for a given hybrid sample of interest maps alleles inherited from potential parents
### Info:
perl VCF-Public-Private-Alleles-Script.pl

Usage: perl VCF-Public-Private-Alleles-Script.pl -t <name_tab_lineage_tab_location> -v <vcf_file_individual_of_interest> -o <output_directory_name>

--------------------------------------------------------

Script to find parental lineages of hybrid isolates

--------------------------------------------------------

Example of tab file for -t:
    04CN-64-029  VNI 04CN-64-029.vcf
E.g.,

perl VCF-Public-Private-Alleles-Script.pl -t name_location.tab -v 881205_filtered.vcf -o test_for_script


## Mapping_Allele_Inheritance_10kb_Windows_Script.pl
This script uses the output chromosome and position tab outputs from VCF-Public-Private-Alleles-Script.pl. The outputs are the tab files which are containing allele inheritance information for the likeli parent
### Info:
Usage:
    perl Mapping_Allele_Inheritance_10kb_Windows_Script.pl -r <reference fasta file> -f <first parent tab file> -s <second parent tab file> -o <output file prefix>

Example:
    perl Mapping_Allele_Inheritance_10kb_Windows_Script.pl -r reference.fa -f parent1.tab -s parent2.tab -o output_prefix

Options:
    -r <reference fasta file>    : Input reference fasta file
    -f <first parent tab file>   : Input tab-delimited file for first parent
    -s <second parent tab file>  : Input tab-delimited file for second parent
    -o <output file prefix>      : Prefix for output files




