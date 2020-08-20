#!/usr/bin/env perl

use warnings;
use strict;

use Carp;
use Data::Dumper;

sub usage {
    print "Syntax: $0 <vcf.map> <sumstats>\n";
    print "Checks and updates a sumstats file with RefPanelAF, INFO and Genotyped columns. Writes to stdout\n";
    exit 1;
}

sub cmp_bp {
    my $left = shift;
    my $right = shift;
    my @l = split(/:/, $left);
    my @r = split(/:/, $right);

    return -1 if($l[0] < $r[0]);
    return 1 if($l[0] > $r[0]);
    return -1 if($l[1] < $r[1]);
    return 1 if($l[1] > $r[1]);
    return 0;
}

usage unless @ARGV == 2;

my $vcfmap_file = shift @ARGV;
my $sumstats_file = shift @ARGV;

# Parse vcf.map
# 1 16833070 rs696608 A G . PASS RefPanelAF=0.590222;AN=16550;AC=14194;INFO=0.725983


# Read sumstats
my @sumstats;
open my $sumstats_fh, '<', $sumstats_file or die($!);
my $header = <$sumstats_fh>;
while(<$sumstats_fh>) {
    chomp;
    my @parts = split(/\s+/);
    my $partsref = \@parts;
    push(@sumstats, $partsref);
}


# Read vcf.map
my %vcfentries;
open my $vcfmap_fh, '<', $vcfmap_file or die($!);
while(<$vcfmap_fh>) {
    chomp;
    my @parts = split(/\s+/);
    my $info = $parts[7];

    my $typed = 0;
    my $af = 0.0;
    my $infoscore = 0.0;
    my $refa1 = '';

    $refa1 = $parts[3];
    $typed = 1 if $info =~ /(TYPED|PHASED)/;
    $af = $2 if $info =~ /(RefPanelAF|RefPAF)=(\d+\.\d+);/;
    $infoscore = $2 if $info =~ /(R2|INFO)=(\d+(\.\d+)?)/;

    my $vcfkey = $parts[0].":".$parts[1]."_".$parts[3].$parts[4];

    #    print "$vcfkey\n";
    $vcfentries{$vcfkey} = [ $typed, $af, $infoscore, $refa1 ];
}

sub cmp_sumstats {
    my @pa = @$a;
    my @pb = @$b;

    return $pa[0] <=> $pb[0] unless $pa[0] == $pb[0];
    return $pa[1] <=> $pb[1];
}

chomp $header;
$header .= "\t" . join("\t", qw(R2 RefPanelAF_A1 RefPanelAF Genotyped));
print "$header\n";
my @sorted = sort cmp_sumstats @sumstats;

foreach(@sorted) {
    my @orig_parts =@$_;
    print join("\t", @orig_parts);
    my $key = "$orig_parts[0]:$orig_parts[1]_$orig_parts[4]$orig_parts[3]";
    my $altkey = "$orig_parts[0]:$orig_parts[1]_$orig_parts[3]$orig_parts[4]";
    #    print $key . "\n";
    my $partref = $vcfentries{$key};
    $partref = $vcfentries{$altkey} unless defined($partref);

    if(defined($partref)) {
        my @newparts = @{$partref};
        print "\t". join("\t",($newparts[2], $newparts[3],  $newparts[1])) . "\t";
        if($newparts[0]) {
            print "GenotypedImputed";
        } else {
            print "ImputedOnly";
        }
    } else {
        print "\t". join("\t", qw(NA NA NA GenotypedOnly));
    }

    print "\n";
}
