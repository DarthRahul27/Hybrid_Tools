#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Std;


# Command line options
our ($opt_r, $opt_f, $opt_s, $opt_o);
getopts('r:f:s:o:');

# Usage information
my $usage = <<END_USAGE;
Usage:
    perl $0 -r <reference fasta file> -f <first parent tab file> -s <second parent tab file> -o <output file prefix>

Example:
    perl $0 -r reference.fa -f parent1.tab -s parent2.tab -o output_prefix

Options:
    -r <reference fasta file>    : Input reference fasta file
    -f <first parent tab file>   : Input tab-delimited file for first parent
    -s <second parent tab file>  : Input tab-delimited file for second parent
    -o <output file prefix>      : Prefix for output files
END_USAGE

# Check for required options
die $usage unless ($opt_r && $opt_f && $opt_s && $opt_o);

# Assign command line arguments to variables
my $reference_file     = $opt_r;
my $first_parent_file  = $opt_f;
my $second_parent_file = $opt_s;
my $output_prefix      = $opt_o;

# Step 1: Generate 10kb windows from reference fasta file
my @windows = making_10kb_windows($reference_file);

# Step 2: Read tab files and tally positions across genome
my $first_parent_positions  = read_tab_file($first_parent_file);
my $second_parent_positions = read_tab_file($second_parent_file);

# Step 3: Calculate tally of positions in each window for both parents
my %tally_data = calculate_window_tally(\@windows, $first_parent_positions, $second_parent_positions);

# Step 4: Output data to tab-delimited files for R plotting
my $output_file = output_data_to_file(\%tally_data, $output_prefix);


print "Aggregated data written to $output_prefix.all_chromosomes.position_tally.tab\n";

# Subroutine to generate 10kb windows from reference fasta file
sub making_10kb_windows {
    my ($reference_file) = @_;
    my $window_size = 10000;
    my @windows;

    my $seq_in = Bio::SeqIO->new(-file => $reference_file, -format => 'fasta');
    while (my $seq = $seq_in->next_seq) {
        my $seq_id = $seq->id;
        my $seq_length = $seq->length;

        for (my $start = 0; $start < $seq_length; $start += $window_size) {
            my $end = $start + $window_size - 1;
            $end = $seq_length - 1 if $end >= $seq_length;
            my @window_data = ($seq_id, $start, $end);
            push @windows, \@window_data;
        }
    }

    return @windows;
}

# Subroutine to read tab-delimited file and tally positions within a window
sub read_tab_file {
    my ($file) = @_;
    my %positions;

    open my $fh, '<', $file or die "Cannot open $file: $!";
    <$fh>; # Skip header

    while (my $line = <$fh>) {
        chomp $line;
        my ($chromosome, $position) = split /\t/, $line;
        push @{$positions{$chromosome}}, $position;
    }

    close $fh;
    return \%positions;
}

# Subroutine to calculate tally of positions within each window for both parents
sub calculate_window_tally {
    my ($windows_ref, $first_parent_positions, $second_parent_positions) = @_;
    my %tally_data;

    foreach my $window_ref (@$windows_ref) {
        my ($seq_id, $start, $end) = @$window_ref;

        my $count_first_parent  = tally_positions_in_window($first_parent_positions->{$seq_id}, $start, $end);
        my $count_second_parent = tally_positions_in_window($second_parent_positions->{$seq_id}, $start, $end);

        push @{$tally_data{$seq_id}}, {
            Start              => $start,
            End                => $end,
            First_Parent_Count  => $count_first_parent,
            Second_Parent_Count => $count_second_parent,
            Running_Position   => ($start + $end) / 2
        };
    }

    return %tally_data;
}

# Subroutine to tally positions within a specified window
sub tally_positions_in_window {
    my ($positions, $start, $end) = @_;
    my $count = 0;

    foreach my $pos (@$positions) {
        if ($pos >= $start && $pos <= $end) {
            $count++;
        }
    }

    return $count;
}

# Subroutine to output data to tab-delimited file for R plotting
sub output_data_to_file {
    my ($tally_data_ref, $output_prefix) = @_;

    my $all_chromosomes_tab = "$output_prefix.all_chromosomes.position_tally.tab";
    open my $all_tab_fh, '>', $all_chromosomes_tab or die "Cannot open $all_chromosomes_tab for writing: $!";
    print $all_tab_fh "Chromosome\tStart\tEnd\tRunning_Position\tFirst_Parent_Count\tSecond_Parent_Count\n";

    foreach my $seq_id (keys %$tally_data_ref) {
        foreach my $data (@{$tally_data_ref->{$seq_id}}) {
            print $all_tab_fh "$seq_id\t$data->{Start}\t$data->{End}\t$data->{Running_Position}\t$data->{First_Parent_Count}\t$data->{Second_Parent_Count}\n";
        }
    }
    close $all_tab_fh;
}


