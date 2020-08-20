#!/usr/bin/env perl

use warnings;
use strict;

use 5.020;

use Data::Dumper;

die("Syntax: $0 <file.bim> <file.annotations>") unless @ARGV == 2;

my $bimfile = $ARGV[0];
my $annotfile = $ARGV[1];

my %annotations;

# parse annotations
# ichip_id -> (rs_id, alleleA, alleleB, strand, pos_hg18, pos_hg19)

open(my $annot_fh, "<$annotfile") or die("Could not open $annotfile: $!");
while(<$annot_fh>) {
    chomp;
    my @l = split(' ');
    $annotations{$l[0]} = { 'rs' => $l[1], 'a1' => $l[2], 'a2' => $l[3], 'strand' => $l[4], 'hg18' => $l[5], 'hg19' => $l[6] };
}
close $annot_fh;

my %complement = (
    A => 'T',
    T => 'A',
    C => 'G',
    G => 'C',
    D => 'D',
    I => 'I'
    );


open(my $bim_fh, "<$bimfile") or die("Could not open $bimfile: $!");
while(<$bim_fh>) {
    chomp;
    my @l = split(/\s+/);
    my $ichip = $l[1];
#    say "Looking up $ichip";
    my %ann = %{$annotations{$ichip}};

    my $l4 = $l[4];
    my $l5 = $l[5];

    # was: checkalleles (fh, ichip, l4, l5, a1, a2, strand)
    if($ann{strand} eq '-') {
        $l4 = $complement{$l4};
        $l5 = $complement{$l5};
    }
#    say "Alleles: $l4 $l5 $ann{a1} $ann{a2}";
    next if $ann{strand} eq '---';

    next if  (($l4 eq $ann{a1} or $l4 eq '0') and ($l5 eq $ann{a2} or $l5 eq '0')) or (($l4 eq $ann{a2} or $l4 eq '0') and ($l5 eq $ann{a1} or $l5 eq '0'));


    if ( (($l4 eq '0' or $complement{$l4} eq $ann{a1}) and ($l5 eq '0' or $complement{$l5} eq $ann{a2})) or (($l4 eq '0' or $complement{$l4} eq $ann{a2}) and ($l5 eq '0' or $complement{$l5} eq $ann{a1})) ) {
        say $ichip;
    } else {
        say STDERR "$ichip has wrong alleles in bim file";
    }
}
close $bim_fh;

