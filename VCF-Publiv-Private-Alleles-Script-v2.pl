#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;
use File::Path 'make_path';

## RA + RF Script to find parental lineages 

my $usage = "Usage: perl $0 -t <name_tab_lineage_tab_location> -v <vcf_file_individual_of_interest> -o <output_directory_name>

--------------------------------------------------------

Script to find parental lineages of hybrid isolates

--------------------------------------------------------

Example of tab file for -t:
    04CN-64-029  VNI 04CN-64-029.vcf
    ERS1238896   VNB ERS1238896.vcf
";

# Parse command line arguments
our ($opt_t, $opt_v, $opt_o);
getopts('t:v:o:');
die $usage unless ($opt_t && $opt_v && $opt_o);

my $output_directory = $opt_o;
my $tab_file = $opt_t;
my $vcf_file_of_interest = $opt_v;

# Extract the sample name from the VCF file name
my ($sample_of_interest) = $vcf_file_of_interest =~ /([^\/]+)\.vcf$/;

# Create output directory if it doesn't exist
unless (-d $output_directory) {
    make_path($output_directory) or die "Failed to create directory $output_directory: $!";
}

# Parse tab file
my %file_data_info = parse_tab_file($tab_file);

# Store SNP information for each lineage
my %lineage_snp_info;
foreach my $name (keys %file_data_info) {
    my $lineage = $file_data_info{$name}{lineage};
    my $vcf_file = $file_data_info{$name}{location};
    my %snp_info = save_snp_info($vcf_file);
    $lineage_snp_info{$lineage} = \%snp_info;
}

# Find and save public allele information
find_public_allele_info($vcf_file_of_interest, $output_directory, \%lineage_snp_info, $sample_of_interest);

## Subroutine to parse tab file
sub parse_tab_file {
    my ($tab_file) = @_;
    my %file_data;

    warn "Parsing tab file...";

    open my $tfh, '<', $tab_file or die "Cannot open tab file $tab_file: $!";
    while (my $line = <$tfh>) {
        chomp $line;

        my ($name, $lineage, $location) = split /\t/, $line;

        if (!defined $lineage || !defined $location || !-e $location) {
            warn "Information for sample $name is not complete or location file $location does not exist";
            die;
        } else {
            warn "$name info: lineage $lineage and location $location exists";
        }

        $file_data{$name} = {lineage => $lineage, location => $location};
    }

    close $tfh;
    warn "Parsing tab file complete";
    return %file_data;
}

## Subroutine to save SNP info from VCF file
sub save_snp_info {
    my ($vcf_file) = @_;
    my %snp_info;

    warn "Saving SNP info for VCF file: $vcf_file";

    open my $vcf_fh, '<', $vcf_file or die "Cannot open VCF file $vcf_file: $!";

    while (my $line = <$vcf_fh>) {
        chomp $line;
        next if $line =~ /^#/;

        my @fields = split /\t/, $line;
        my ($chr, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @sample) = @fields;

        my ($genotype) = (split /:/, $sample[0]);

        next if $genotype eq '0';
        next if length($ref) != 1 or length($alt) != 1;

        if ($genotype eq '1') {
            push @{$snp_info{$chr}{$pos}}, $alt;
        }
    }

    close $vcf_fh;
    warn "SNP info saved for VCF file: $vcf_file";
    return %snp_info;
}

## Subroutine to find public allele information
sub find_public_allele_info {
    my ($vcf_file, $output_directory, $lineage_snp_info, $sample_of_interest) = @_;
    my %snp_info_for_sample_of_interest = save_snp_info($vcf_file);
    my %snp_isolates;

    warn "Saving public allele information";

    foreach my $lineage (keys %$lineage_snp_info) {
        foreach my $chr (keys %{$lineage_snp_info->{$lineage}}) {
            foreach my $pos (keys %{$lineage_snp_info->{$lineage}{$chr}}) {
                foreach my $alt (@{$lineage_snp_info->{$lineage}{$chr}{$pos}}) {
                    my $snp_key = "$chr:$pos:$alt";
                    $snp_isolates{$snp_key}{$lineage} = 1;
                }
            }
        }
    }

    foreach my $chr (keys %snp_info_for_sample_of_interest) {
        foreach my $pos (keys %{$snp_info_for_sample_of_interest{$chr}}) {
            foreach my $alt (@{$snp_info_for_sample_of_interest{$chr}{$pos}}) {
                my $snp_key = "$chr:$pos:$alt";
                $snp_isolates{$snp_key}{'sample_of_interest'} = 1;
            }
        }
    }

    my $output_file_prefix = "VCF-Allele-Inheritance";
    my %public_alleles;
    my %summary_counts;

    foreach my $snp_key (keys %snp_isolates) {
        next unless exists $snp_isolates{$snp_key}{'sample_of_interest'};

        my @lineage_names = grep { $_ ne 'sample_of_interest' } keys %{$snp_isolates{$snp_key}};
        next unless @lineage_names;

        my $lineage_names_str = join('_', sort @lineage_names);
        my ($chromosome, $position, $alternative) = split /:/, $snp_key;

        my $outfile = "$output_directory/${output_file_prefix}_public_${lineage_names_str}_$sample_of_interest.tab";

        open my $ofh, '>>', $outfile or die "Cannot open $outfile: $!";
        warn "Writing to outfile: $outfile";
        print $ofh "$chromosome\t$position\n";
        close $ofh;

        $summary_counts{$lineage_names_str}++;
    }

    my $summary_file = "$output_directory/${output_file_prefix}_summary_file_$sample_of_interest.tab";
    open my $summary_fh, '>', $summary_file or die "Cannot open $summary_file: $!";
    warn "Writing summary file";
    print $summary_fh "Names\tCounts\n";
    foreach my $lineage_names_str (keys %summary_counts) {
        my $count = $summary_counts{$lineage_names_str};
        print $summary_fh "$lineage_names_str\t$count\n";
    }

    close $summary_fh;
}
