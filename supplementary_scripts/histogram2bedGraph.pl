#!/usr/bin/perl
use strict;

# This script convert a histogram file with the format specified below 
# to a bedGraph file using as score the "counts" column and as coordinates the 
# IS_coord that corresponds to the insertion site.

# The programs also normalize the scores (peaks) as peaks per million peaks

##Chromosome	flank_5	flank_3	counts	IS_coord	read_strand
#chr22	24310207	24310271	1	24310241	+
#chr12	68572516	68572580	1	68572546	-
#chr10	34723937	34724001	8	34723967	-
#chr2	190648002	190648066	1	190648036	+
my %peaks;
my $total_peaks = 0;

while(<>){
    next if m/^#/;
    chomp;
    my ($chr, $end5, $end3, $counts, $is_coord, $strand) = split /\t/;
    $total_peaks += $counts;
    if ($strand eq "+"){
        my $new5 = $is_coord - 1;
        my $name = "$chr:$end5-$end3";
        $peaks{$name} = "$chr\t$new5\t$is_coord\t$counts";
    } else {
        my $new3 = $is_coord + 1;
        my $name = "$chr:$end5-$end3";
        $peaks{$name} = "$chr\t$is_coord\t$new3\t$counts";
    }
}

# Print out normalized peaks
foreach my $name (keys %peaks){
    my ($chr, $end5, $end3, $counts) = split("\t", $peaks{$name});
    my $norm_counts = ($counts / $total_peaks) * 1000000;
    print "$chr\t$end5\t$end3\t$norm_counts\n";
}


