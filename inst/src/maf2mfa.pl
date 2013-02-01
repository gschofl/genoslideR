#!/usr/bin/perl 

# Title: maf2mfa.pl
# Author: Jochen Bick

use strict;
use Data::Dumper;
#use Parallel;
my $seq = $ARGV[0];
my $maln_seq = "";

my $allfile = readall($seq);
my @sequence = split(/\n\n/, $allfile);
my @seq_names = findgenomes($allfile);

#print Dumper(@sequence);
my %hash = hashfiller(@seq_names);
my $count = 1;

foreach my $aln_frag ( @sequence ){
    my @elements = split(/\ns /, $aln_frag);
    shift(@elements);
    my $aln_size;
    foreach my $k ( @elements ){
	my $genome;
	my $alignment;
	#$k =~ s/.+\|(.+)\|\t(\d+)\t(\d+)\t([\+\-])\t\d+\t([\w\-]+)/$1:$2:$3:$4,$5/;
	if($k =~ /.+\|(.+)\|\t(\d+)\t(\d+)\t([\+\-])\t\d+\t([\w\-]+)/){
	    $genome = $1;
	    $alignment = "$1:$2:$3:$4,$5";
	    $aln_size = length($5);
	    
	}
	#print $genome."\n";
	push(@{$hash{$genome}}, $alignment);
	#print $k."\n";
    }
    foreach my $ele (keys %hash){
	#my $nums = @{$hash{$ele}};
	if($hash{$ele}){
	    my $a = $hash{$ele};
	    my $size = @$a;
	    if($size =~ /$count/){
		#print $hash{$ele};
	    }else{
		my $gaps = "-"x$aln_size;
		#print Dumper($gaps);
		push(@{$hash{$ele}}, ",$gaps");
	    }
	}else{
	    my $gaps = "-"x$aln_size;
	    push(@{$hash{$ele}}, ",".$gaps);
	}
    }
    $count++;
}
#print Dumper(%hash);
printhash(\%hash);


sub printhash{
    my %h = %{$_[0]};
   # print Dumper(%h);
    foreach my $ele (sort keys %h){
	my $a = $h{$ele};
#	print @$a;
#	exit;
	my @header;
	my @aln;
	foreach my $e (@$a){
	    my ($c , $d) = split(/,/,$e);
	    push(@header, $c);
	    push(@aln, $d);
	}
	my $printheader = join(" ", @header)."\n";
	$printheader =~ s/  +/ /g;
	print ">$ele".$printheader;
	print join("", @aln)."\n";
    } 
}


sub findgenomes{
    my $file = $_[0];
    my @returnarr;
    $file =~ s/a .*\n//g;
    $file =~ s/s .*\|([\w\d\.]+)\|.+\n*/$1,/g;
    @returnarr = split(/,/, $file);
    @returnarr = keys %{{ map { $_ => 1 } @returnarr }}; # unique array
    return @returnarr;
}


sub hashfiller{
    my @g = @_;
    my %h;
    foreach my $i (@g){
#	$h{$i} = "";
	push(@{$h{$i}}, "");
    }
    return %h;
}


sub readall{
    open(FILE, "<$_[0]") or die "Cannot open $_[0]: $!";
    undef $/;# normaly $/ is set on newline charater
    my $all = (<FILE>);# so now it reads the whole file in one scaler
    unless( length($all) > 0 ) {
	die "Zero length input file\n";
    }
    close(FILE);
    return $all;
}



