#!/usr/bin/env perl

use warnings;
use strict;

use Data::Dumper;

use 5.020;

sub  trim {
    my $s = shift;
    $s =~ s/^\s+|\s+$//g;
    return $s;
};

my %sections;
my $section;

sub usage {
    say "Usage: config-convert.pl [input-config] [export-section]";
    say "\tExample: config-convert.pl ../PoAtools_2/config mod_rs";
    say "\twill export the section named 'mod_rs' to stdout.";
    exit;
}

usage unless @ARGV == 2;

open INFILE, '<', $ARGV[0] or die($!);
while(<INFILE>) {
    chomp;

    next if /^#/;
    next if /^$/;

    if(/^\[(.*)\]/) {
        $section = $1;
#        say $1;
        $sections{$section} = {};
        next;
    }

    my @parts = split(':');
    if(@parts == 1) {
        @parts = split('=');
    }
    $sections{$section}{trim($parts[0])} = trim($parts[1]);
}

say "// -*- mode:groovy -*-";
say "// Exporting section $ARGV[1]...";
say "params {";
my ($k, $v);
while(my ($k, $v) = each %{$sections{$ARGV[1]}}) {
    say "  $k = \"$v\"";
}
say "}";
