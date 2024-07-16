# Hybrid_Isolate_Allele_Inheritance_Script

## This script uses VCF files to get single nucleotide polymorphism data, assigning those snps to designated lineage and for a given hybrid sample of interest maps alleles inherited from potential parents

### Main Script:
perl VCF-Publiv-Private-Alleles-Script-v2.pl 

Usage: perl VCF-Publiv-Private-Alleles-Script-v2.pl -t <name_tab_lineage_tab_location> -v <vcf_file_individual_of_interest> -o <output_directory_name>

--------------------------------------------------------

Script to find parental lineages of hybrid isolates

--------------------------------------------------------

Example of tab file for -t:
    04CN-64-029  VNI 04CN-64-029.vcf
E.g.,

perl VCF-Publiv-Private-Alleles-Script-v2.pl -t name_location.tab -v 881205_filtered.vcf -o test_for_script
