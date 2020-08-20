#!/usr/bin/env perl

use warnings;
use strict;

use 5.020;

use Data::Dumper;
use Switch 'fallthrough';

die("Syntax: $0 <format> <file.bim> <file.annotations> <from> <to>") unless @ARGV == 5;

my $format = $ARGV[0];
my $bimfile = $ARGV[1];
my $annotfile = $ARGV[2];
my $chip_build = $ARGV[3];
my $switch_to = $ARGV[4];

my %annotations;

# parse annotations
# ichip_id -> (rs_id, alleleA, alleleB, strand, pos_hg18, pos_hg19)

open(my $annot_fh, "<$annotfile") or die("Could not open $annotfile: $!");
while(<$annot_fh>) {
    chomp;
    my @l = split(' ');

#    switch($format) {
#        case 'Immunochip' {
            $annotations{$l[0]} = { 'rs' => $l[1], 'a1' => $l[2], 'a2' => $l[3], 'strand' => $l[4], 'hg18' => $l[5], 'hg19' => $l[6] };
#            #$annotations{$l[3]} = { 'rs' => $l[4], 'a1' => $l[7], 'a2' => $l[8], 'strand' => $l[5], 'hg18' => $l[1], 'hg19' => $l[2] };#
#
 #           last;
  #      }
#        case 'Exomechipv1' {}
#        case 'Exomechipc1-1' {}
#        case 'HumanCoreExome24v1' {}
#        case 'GSAarrayv1' {
#            $annotations{$l[2]} = { 'rs' => $l[3], 'a1' => $l[6], 'a2' => $l[7], 'strand' => $l[4], 'hg18' => $l[1], 'hg19' => $l[1] };
#            last;
#        }
#    }
}
close $annot_fh;


open(my $bim_fh, "<$bimfile") or die("Could not open $bimfile: $!");
while(<$bim_fh>) {
    chomp;
    my @l = split(/\s+/);
    my $ichip = $l[1];

    if(defined $annotations{$ichip}) {
        my %ann = %{$annotations{$ichip}};

        if ($ann{hg19} eq 'NA') {
            $ann{hg19} = '-9999999';
        }

        if ($chip_build eq 'hg18' and $switch_to eq 'hg19') {
            say "$l[0]\t$ann{rs}\t$l[2]\t$ann{hg19}\t$l[4]\t$l[5]";
        } else {
            say "$l[0]\t$ann{rs}\t$l[2]\t$l[3]\t$l[4]\t$l[5]";
        }
    } else {
        say STDERR "Could not find SNP $ichip in annotation file!";
    }

}
close $bim_fh;

