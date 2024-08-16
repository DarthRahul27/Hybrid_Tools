#!/usr/bin/perl -w 
use strict;
use Getopt::Std;
use Data::Dumper;
use File::Path 'make_path';

my $usage = "Usage: perl $0 -t <name_tab_lineage_file> -v <vcf_file_list> -o <outptut_directory_prefix>


Script to find parental lineages of hybrid isolates

--------------------------------------------------------

Example of tab file for -t:
    04CN-64-029  VNI 04CN-64-029.vcf
    ERS1238896   VNB ERS1238896.vcf

Example of VCF list file for -v:
    vcf_file1.vcf
    vcf_file2.vcf
    ...

";

#Parsing command line arguements 
our ($opt_t, $opt_v, $opt_o);
getopts('t:v:o:');
die $usage unless ($opt_t && $opt_v && $opt_o);
my $output_directory = $opt_o;
my $tab_file = $opt_t;
my $vcf_file_list = $opt_v;

#Create output directory if it doesnt exist
unless (-d $output_directory) {
	make_path($output_directory) or die "Cannot make $output_directory: $!";
}

#Parsing tab file 
my $file_data_info = parse_tab_file($tab_file);

my %lineage_snp_info;
foreach my $name (keys %{$file_data_info}) {
	my $lineage = $$file_data_info{$name}{'lineage'};
	my $vcf_file = $$file_data_info{$name}{'location'};
	my $snp_info = save_snp_info($vcf_file);
	$lineage_snp_info{$lineage} = $snp_info;
}


my @vcf_files = parse_vcf_file_list($vcf_file_list);


hybrid_or_not(\@vcf_files, \%lineage_snp_info, $output_directory);

warn "Script Complete";




## Subroutine to parse tab file
sub parse_tab_file {
    my ($tab_file) = @_;

    warn "parse_tab_file: $tab_file\n";

    my %file_data;
    open my $tfh, '<', $tab_file or die "Cannot open tab file $tab_file: $!";
    while (my $line = <$tfh>) {
        chomp $line;

        my ($name, $lineage, $location) = split /\t/, $line;

        if (!defined $lineage || !defined $location || !-e $location) {
            die "Information for sample $name is not complete or location file $location does not exist";
        } else {
            warn "$name info: lineage $lineage and location $location exists";
        }
        $file_data{$name}{'lineage'} = $lineage;
        $file_data{$name}{'location'} = $location;
    }
    close $tfh;
    warn "Parsing tab file complete\n";
    return \%file_data;
}

sub save_snp_info {
    my ($vcf_file) = @_;

    warn "save_snp_info: $vcf_file\n";

    my %snp_info;
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
    warn "save_snp_info: finished saving $vcf_file\n";
    return \%snp_info;
}



sub parse_vcf_file_list {
	my ($vcf_file_list) = @_;

	warn "parse_vcf_file_list $vcf_file_list...";

	open my $vcf_list_fh, '<', $vcf_file_list or die "cannot open $vcf_file_list";
	my @vcf_files_of_interest;

	while (my $line = <$vcf_list_fh>) {
		chomp $line;
		push @vcf_files_of_interest, $line;
	}

	close $vcf_list_fh;

	warn "parse_vcf_file_list completed for $vcf_file_list";

	return @vcf_files_of_interest;

}


sub hybrid_or_not {
    my ($vcf_files_ref, $lineage_snp_info_ref, $output_directory) = @_;
    my @vcf_files = @$vcf_files_ref;

    my $summary_file = "$output_directory/VCF-Allele-Inheritance_summary.tab";
    open my $summary_fh, '>', $summary_file or die "Cannot open $summary_file: $!";
    print $summary_fh "VCF File\tTop Lineage\tShared Alleles\tTotal Alleles\tShared Percentage\n";

    foreach my $vcf_file (@vcf_files) {
        my $snp_info_for_sample_of_interest = save_snp_info($vcf_file);

        # Extract the sample name from the VCF file name
        my ($sample_of_interest) = $vcf_file =~ /([^\/]+)\.vcf$/;

        my %snp_isolates;
        foreach my $lineage (sort keys %$lineage_snp_info_ref) {
            foreach my $chr (sort keys %{$$lineage_snp_info_ref{$lineage}}) {
                foreach my $pos (sort { $a <=> $b } keys %{$$lineage_snp_info_ref{$lineage}{$chr}}) {
                    foreach my $alt (@{$$lineage_snp_info_ref{$lineage}{$chr}{$pos}}) {
                        my $snp_key = "$chr:$pos:$alt";
                        $snp_isolates{$snp_key}{$lineage} = 1;
                    }
                }
            }
        }

        foreach my $chr (keys %$snp_info_for_sample_of_interest) {
            foreach my $pos (keys %{$$snp_info_for_sample_of_interest{$chr}}) {
                foreach my $alt (@{$$snp_info_for_sample_of_interest{$chr}{$pos}}) {
                    my $snp_key = "$chr:$pos:$alt";
                    $snp_isolates{$snp_key}{'sample_of_interest'} = 1;
                }
            }
        }

        my %lineage_shared_count;
        my $total_shared_count = 0;
        foreach my $snp_key (keys %snp_isolates) {
            next unless exists $snp_isolates{$snp_key}{'sample_of_interest'};

            my @lineage_names = grep { $_ ne 'sample_of_interest' } keys %{$snp_isolates{$snp_key}};
            next unless @lineage_names == 1;  # Only consider SNPs shared with exactly one lineage

            my $lineage = $lineage_names[0];
            $lineage_shared_count{$lineage}++;
            $total_shared_count++;
        }

        # Find the lineage with the highest shared allele count
        my ($top_lineage) = sort { $lineage_shared_count{$b} <=> $lineage_shared_count{$a} } keys %lineage_shared_count;
        my $top_lineage_count = $lineage_shared_count{$top_lineage} // 0;

        # Calculate the shared percentage based on the total shared allele count
        my $shared_percentage = $total_shared_count > 0 ? ($top_lineage_count / $total_shared_count) * 100 : 0;

        # Output the detailed SNP information
        my $output_file_prefix = "VCF-Allele-Inheritance";
        my $outfile = "$output_directory/${output_file_prefix}_public_${top_lineage}_$sample_of_interest.tab";
        open my $ofh, '>', $outfile or die "Cannot open $outfile: $!";
        foreach my $snp_key (keys %snp_isolates) {
            next unless exists $snp_isolates{$snp_key}{'sample_of_interest'};
            next unless exists $snp_isolates{$snp_key}{$top_lineage};

            my ($chromosome, $position, $alternative) = split /:/, $snp_key;
            print $ofh "$chromosome\t$position\n";
        }
        close $ofh;

        # Write summary file
        print $summary_fh "$vcf_file\t$top_lineage\t$top_lineage_count\t$total_shared_count\t$shared_percentage\n";
    }

    close $summary_fh;
}









