#!/Users/lorenziha/miniconda3/envs/TK_58/bin/perl
use strict;

my $usage = "$0 -f <fasta file> -c <chromosome id: e.g chr2-chr5-chr8 [def=all]>\n\n";
my %arg = @ARGV;
die $usage unless $arg{-f} ;
my $CHR_ID = $arg{-c} || 'all';

my (%h, $id);
open (FASTA, "<$arg{-f}") || die "ERROR, I cannot open $arg{-f}: $!\n\n";
while(<FASTA>){
	chomp;
	if(m/^>(\S+)/){
		$id=$1; 
	}
	elsif(m/^[ACGTNXacgtnx]+$/){
		$h{$id}.= $_;
	}
}; 
close FASTA;

# Search TTAA motif and its position on chromosome
# and print a bed file with the choords

my @x;
if ($CHR_ID eq 'all'){
	@x = keys %h;
} else {
	@x = split(/-/,$CHR_ID);
}
foreach my $id (@x){
	print STDERR "Processing $id\n";
	while($h{$id} =~ m/(TTAA|AATT)/g){
		my $end5 = length($`);
		my $end3 = $end5 + 1;
		my $posid = "$id:$end5-$end3";
		print "$id\t$end5\t$end3\t$posid\n";
	}
	#print STDERR "Processing $id reverse\n";
	#while($h{$id} =~m/(AATT)/g){
        #        my $end5 = length($`);
        #        my $end3 = $end5 + 1;
        #        my $posid = "$id:$end5-$end3";
        #        print "$id\t$end5\t$end3\t-\t$posid\n";
        #}
}





