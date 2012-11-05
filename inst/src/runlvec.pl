#!/usr/bin/perl -w

# Title: rnazpaser.pl
# Author: Jochen Bick

use strict;

my $seq = $ARGV[0];
my $header = $ARGV[1];
#my $maln_seq = <STDIN>;
#my $w_start = $ARGV[1];
my $maln_seq;
#my $maln_seq = "ATGAGCTCATTAGCGTTGAGTCGACGACCTCGTAGAAATAGAAGAACTGAGGCAATTCGTGATTTGGTATCTGAAACTTCTTTATTACCTCAGGATTTCATCTGTCCTTTTTTTGTTAAGGAAGGAAAAAACATACGTGAGGAAATAGAAAGTCTTACAGGTGTGTATAGATGGAGTATAGATCTTCTTTTAAAAGAGATCGAAAGGTTGTGTTCTTTGGGGTTAAGAGCTGTGATTCTTTTTCCCGTCATTCCGAGTCATCTTAAAGATGCTTACGGTTCATACTCTTCTAATCCTAAAAATATTCTATGTAAGAGCATCTATGAAGTAAAAAAAGCTTTTCCTAACTTGTGTGTAATCAGCGACATTGCTTTAGATCCTTATACCACTCATGGTCATGACGGTATTATTG-ATC--GTGGGGAAGTATTGAATGATGAGAGCGTGCGTATATTCGGTAACATAGCTACGTTGCATGCTGAGATGGGCGCTGATGTTGTAGCTCCGAGTGATATGATGGACGGCAGGGTAGCTCATATCCGCTCGAAGTTAGATCAAGCAGGATGGACTCAAACCTTAATTCTCTCTTATAGCGTTAAGTATGCTTCTGCTTTGTATAATCCTTTTCGAGATGCTTTAGGATCTCATTTGCAATCAGGAGATAAACGCAATTATCAGATGAATCCTAAAAATGTTTTAGAAGCTTTACTAGAGTGTTCTTTAGATGAGCAAGAGGGCGCTGATATGCTCATGATAAAACCCGCAGGATTATATCTTGATGTGTTGCATCGAGTGAAAAATAGTACAACATTACCCTTAGCAGCCTATCAAGTCAGTGGTGAATATGCTATGATAGCAGCAGCTTCTACTATGGGGTGGTTGGATAGAGAAAAGATAGTGTATGAATCTTTGATAGCGATAAAACGCGCAGGAGCCGACATGATCATTTCTTACGCAACTCCGTTAATTTTAGAAATGATAGCTTCT---------------TCGAG------AGTGTAA";


my $count = 0;
my $read_start = 0;
open(SEQUENCE, "<$seq") or die "cannot open file: $!\n";
while(<SEQUENCE>){
    if($_ =~ />/){
	$count++;
	if($count == $header){
	    $read_start = 1;
	}elsif($count > $header){
	     last;
	}	    
    }else{
	next if(!$read_start);
	chomp($_);
	$maln_seq .= $_;
    }
}

close(SEQUENCE);

my $gap_pos = 0;
my $gap_length = 0;
my $gap_start = 0;
my $gap_check = 1;
my $start_gap_check = 0;
#print length($maln_seq)."\n";
#open(TEMPOUT, ">/home/jbick/src/R/genoslider2/R/temp.out") or die "cannot open file: $!\n";
my @arr = split(/(?=.)/, $maln_seq);
#my %hash;

foreach my $ele (@arr) {
    $gap_pos++;
    if($ele =~ /\-/){
	if($gap_check){
	    $gap_start = $gap_pos;
	    $gap_check = 0;
	    $start_gap_check = 1;
	}
	$gap_length++;
    }else{
	if($start_gap_check){
	    print $gap_start."\t".$gap_length."\n";
	    $gap_length = 0;
	    $gap_check = 1;
	    $start_gap_check = 0;
	}
	next;
    }
}

if(!$gap_check){
    print $gap_start."\t".$gap_length."\n";
}
 
#print Dumper(%hash);
#print STDOUT join(" ", sort {$a <=> $b} (keys (%hash)));
#print STDOUT join(",", sort {$a <=> $b} (values (%hash)));
#my @values;
#foreach my $i (sort {$a <=> $b} (keys (%hash))){
#    print "$i\t$hash{$i}\n";
   # push(@values, $hash{$i});    
#}
#print "\n";
#print STDOUT join(" ", @values);
#print $gen_p[149];
#print "0,".join(",",@gen_p).";"."1,".join(",",@maln_p);
#print TEMPOUT join("," ,@gen_p);#.";".join(",",@maln_p);
#close(TEMPOUT);
