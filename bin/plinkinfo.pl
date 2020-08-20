#!/usr/bin/env perl

use strict;
use warnings;

use 5.020;
use Switch;
use Data::Dumper;

sub usage {
    say "Usage: plinkinfo.pl <bim> <fam>";
    say "\tExtracts summary information from a Plink file set.";
    exit;
}

usage() unless @ARGV == 2;

my $arg_bim = shift @ARGV;
my $arg_fam = shift @ARGV;

sub parse_fam {
    my $f = shift;
    my %dat;

    $dat{'male-cases'} = 0;
    $dat{'female-cases'} = 0;
    $dat{'unk-cases'} = 0;
    $dat{'male-controls'} = 0;
    $dat{'female-controls'} = 0;
    $dat{'unk-controls'} = 0;
    $dat{'male-unk'} = 0;
    $dat{'female-unk'} = 0;
    $dat{'unk-unk'} = 0;

    open my $fam_fh, '<', $f or die("Could not open $f: $!");
    while(<$fam_fh>) {
        chomp;
        my @parts = split /\s+/;

        switch($parts[4]) { # sex code
            case 1 {
                switch($parts[5]) {
                    case 1 { $dat{'male-controls'}++; }
                    case 2 { $dat{'male-cases'}++; }
                    else { $dat{'male-unk'}++; }
                }
            }

            case 2 {
                switch($parts[5]) {
                    case 1 { $dat{'female-controls'}++; }
                    case 2 { $dat{'female-cases'}++; }
                    else { $dat{'female-unk'}++; }
                }
            }

            else {
                switch($parts[5]) {
                    case 1 { $dat{'unk-controls'}++; }
                    case 2 { $dat{'unk-cases'}++; }
                    else { $dat{'unk-unk'}++; }
                }
            }
        }
    }

    $dat{'females'} = $dat{'female-cases'} + $dat{'female-controls'} + $dat{'female-unk'};
    $dat{'males'} = $dat{'male-cases'} + $dat{'male-controls'} + $dat{'male-unk'};
    $dat{'unknown-sex'} = $dat{'unk-cases'} + $dat{'unk-controls'} + $dat{'unk-unk'};
    $dat{'cases'} = $dat{'male-cases'} + $dat{'female-cases'} + $dat{'unk-cases'};
    $dat{'controls'} = $dat{'male-controls'} + $dat{'female-controls'} + $dat{'unk-controls'};
    $dat{'unknown-pheno'} = $dat{'male-unk'} + $dat{'female-unk'} + $dat{'unk-unk'};
    $dat{'samples'} = $dat{'cases'} + $dat{'controls'} + $dat{'unknown-pheno'};
    %dat
}


sub parse_bim {
    my $b = shift;
    my $lines = 0;
    open my $bim_fh, '<', $b or die("Could not open $b: $!");
    $lines++ while <$bim_fh>;
    $lines
}

my %d = parse_fam($arg_fam);
$d{'variants'} = parse_bim($arg_bim);

foreach my $key (keys %d) {
    my $v = $d{$key};
    say "$key $v";
}

