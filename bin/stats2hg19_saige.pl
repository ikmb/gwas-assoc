#!/usr/bin/env perl

use warnings;
use strict;

# Process CLI args
my $dict_name = shift;
my $stats_name = shift;

# Read dictionary
my %dict;
print "Reading dictionary...\n";
open my $dict_fh, '<', $dict_name or die($!);
my @parts;
my @varparts;

my $count_changer=0;
my $count_keys = 0;

while(<$dict_fh>) {
    chomp;
    @parts = split ' ';
    @varparts = split ':', $parts[4];
    # Insert only if no chromosomes are changed
    if($parts[0] eq $varparts[0]) {
        $dict{"$parts[0]:$parts[3]"} = $parts[2];
        $count_keys++;
    } else {
        $count_changer++;
    }
}
my $count_dict = keys %dict;
print "  Skipping $count_changer chromosome-changing rules.\n";
print "  Added $count_keys keys to dictionary.\n";
print "  $count_dict unique keys in dictionary.\n";
print "Processing stats...\n";
open my $stats_fh, '<', $stats_name or die($!);
open my $outfh, '>', $stats_name.".b37" or die($!);

my $header;
$header = <$stats_fh>;

print $outfh $header;
my $count_skipped=0;
while(<$stats_fh>) {
    chomp;
    @parts = split ' ';
    if(exists($dict{"$parts[0]:$parts[1]"})) {
        $parts[1] = $dict{"$parts[0]:$parts[1]"};
        print $outfh join(' ', @parts) . "\n";
    } else {
#        print "parts: @parts\n";
        $count_skipped++;
    }
}
print "  Of $. stats results, $count_skipped could not be mapped.\n";

