#!/usr/bin/env perl

use Data::Dumper;

use strict;

sub Usage {
    print "Syntax: $0 <hwe-input> <filter-index> <whole-thresholds> <whole-outliers> <worstbatchremoved-outliers> <perbatch-thresholds> <perbatch-outliers>\n";
    print "\thwe-input:         output from hwe-calculation script\n";
    print "\tfilter-index:      FDR-filtering index (1 is 10^-2, 2 is 10^-3, n is 10^(-n-1))\n";
    print "\twhole-thresholds:  FDR threshold table containing outlier count and corrected p-val\n";
    print "\twhole-outliers:    outliers wrt. filter-index\n";
    print "\tworstbatchremoved-outliers outliers wrt. filter-index with the worst-performing batch removed\n";
    print "\tperbatch-thresholds: FDR threshold table containing variant counts failing at least one or two batches\n";
    print "\tperbatch-outliers: outliers wrt. perbatch-thresholds and filter-index\n";
    exit 1;
}

# FDR thresholds that false discoveries will be estimated for
my @fdr_thresholds = (1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10);

Usage if @ARGV != 7;


# Process CLI args
my ($hwe_fn, $filter_threshold, $thresholds_fn, $whole_fn, $worstbatchremoved_fn, $perbatch_thresholds_fn, $perbatch_fn) = @ARGV;

# Read HWE result file
my @hwe;
my @hwe_worst_batch_removed;
my @hwe_per_batch;

print STDERR "Loading $hwe_fn...\n";

my $batchcount;
open my $hwefile, '<', $hwe_fn or die($!);
while(<$hwefile>) {
    chomp;
    my @line = split /\t/;
    $batchcount = (@line-6)/4;

    # add whole collection (variant, pval)
    push(@hwe, [ $line[1], $line[4] ]);

    # find worst batch
    my $min_this_batch = 1;
    my $min_not_this_batch = 1;
    for(my $b = 0; $b < $batchcount; $b++) {
        my $this_batch = $line[5+$b*2];
        my $not_this_batch = $line[5+$b*2+1];
        if($this_batch < $min_this_batch) {
            $min_this_batch = $this_batch;
            $min_not_this_batch = $not_this_batch;
        }

        #
        # while we're here, also collect batch-only values
        #

        # initialize counters, if necessary
        if(@hwe_per_batch == 0) {
            for(my $i = 0; $i < $batchcount; $i++) {
                $hwe_per_batch[$i] = [];
            }
        }

        # add batch-wise records (variant, pval)
        push(@{$hwe_per_batch[$b]}, [ $line[1], $this_batch] );
    }

    # add variant with pval of that one collection with the
    # lowest-performing batch removed (i.e. with the lowest batch-exclusive pval)
    push(@hwe_worst_batch_removed, [ $line[1], $min_not_this_batch ]);
}

print STDERR "Found " . scalar(@hwe) . " SNPs in " . $batchcount . " batches.\n";
print STDERR "Estimating FDR for whole collection and worst-batch-removed...\n";

#############################################################################
# Actual processing starts here, according to Benjamini-Hochberg procedure

# sort collections by their respective p-value, ascending
@hwe = sort { $a->[1] <=> $b->[1] } @hwe;
@hwe_worst_batch_removed = sort { $a->[1] <=> $b->[1] } @hwe_worst_batch_removed;
foreach(@hwe_per_batch) {
    @$_ = sort { $a->[1] <=> $b->[1] } @$_;
}

# pre-declare
my $rank = 0.0;
my $pthreshold = 0.0;
my %whole_fdr;
my %worst_removed_fdr;

# Correct p-value for multiple testing, record the corrected p-value theshold
# and the (estimated) number of false discoveries among the whole collection
# and for the whole collection without the worst-perfoming batch.
# Both uncorrected p-values are given in the input file.
for(my $j=0; $j <@fdr_thresholds; $j++) {
    my $whole_threshold_hit = 0;
    my $worst_removed_threshold_hit = 0;

    for(my $i=0; $i <@hwe; $i++) {
        $rank = ($i+1) / (@hwe+1);
        $pthreshold = $rank * $fdr_thresholds[$j];

        if(!$whole_threshold_hit && $hwe[$i]->[1] > $pthreshold) {
            $whole_fdr{$j}{'p'} = $hwe[$i-1]->[1];
            $whole_fdr{$j}{'count'} = $i;
            $whole_threshold_hit = 1;
        }
        if(!$worst_removed_threshold_hit && $hwe_worst_batch_removed[$i]->[1] > $pthreshold) {
            $worst_removed_fdr{$j}{'p'} = $hwe_worst_batch_removed[$i-1]->[1];
            $worst_removed_fdr{$j}{'count'} = $i;
            $worst_removed_threshold_hit = 1;
        }
        last if $whole_threshold_hit and $worst_removed_threshold_hit;
    }
}

print STDERR "Estimating FDR per batch...\n";

my @perbatch_fdr;
for(my $b=0; $b<$batchcount; $b++) {
    push(@perbatch_fdr, [ map { [] } @fdr_thresholds]);
}

# ...and now for each batch itself
for(my $b=0; $b < $batchcount; $b++) {
    for(my $f=0; $f<@fdr_thresholds; $f++) {
        for(my $i=0; $i <@hwe; $i++) {
            $rank = ($i+1) / (@hwe+1);
            $pthreshold = $rank * $fdr_thresholds[$f];
            if($hwe_per_batch[$b]->[$i]->[1] > $pthreshold) {
                $perbatch_fdr[$b]->[$f] = $i-1;
                last;
            }
        }
    }
}

#############################################################################
# Dump corrected p-value thresholds and expected false discovery counts

print STDERR "Writing FDR thresholds to $thresholds_fn...\n";
open my $thresfh, '>', $thresholds_fn or die($!);
print $thresfh "FDR\tFail_allbatches\tHWE_pval_allbatches\tFail_worstbatchremoved\tHWE_pval_worstbatchremoved\n";
for(my $j=0; $j < @fdr_thresholds; $j++) {
    print $thresfh sprintf("%e\t%d\t%e\t%d\t%e\n",
        $fdr_thresholds[$j],
        $whole_fdr{$j}{'count'}, $whole_fdr{$j}{'p'},
        $worst_removed_fdr{$j}{'count'}, $worst_removed_fdr{$j}{'p'}
    );
}
close $thresfh;

#
# Dump variants that are rejected with respect to the corrected p-value
#
print STDERR "Writing whole-collection outliers to $whole_fn...\n";
open my $wholefh, '>', $whole_fn or die($!);
# print $wholefh "Variant\tP_HWE\n";
for(my $i=0; $i < $whole_fdr{$filter_threshold}{'count'}; $i++) {
    print $wholefh sprintf("%s\t%e\n", $hwe[$i]->[0], $hwe[$i]->[1]);
}
close $wholefh;

print STDERR "Writing worst-batch-removed outliers to $worstbatchremoved_fn...\n";
open my $wbmfh, '>', $worstbatchremoved_fn or die($!);
# print $wbmfh "Variant\tP_HWE\n";

for(my $i=0; $i < $worst_removed_fdr{$filter_threshold}{'count'}; $i++) {
    print $wbmfh sprintf("%s\t%e\n", $hwe_worst_batch_removed[$i]->[0], $hwe_worst_batch_removed[$i]->[1]);
}
close $wbmfh;

# Determine per-batch outliers that fail the corrected p-value in at
# least one or two batches
#
print STDERR "Writing per-batch outliers to $perbatch_fn...\n";
open my $perbatchfh, '>', $perbatch_fn or die($!);
my %excludes_fail_one;
my %excludes_fail_two;


for(my $b=0; $b < $batchcount; $b++) {
    for(my $i=0; $i < $perbatch_fdr[$b]->[$filter_threshold]; $i++) {
        my $variant = $hwe_per_batch[$b]->[$i]->[0];
        if(not exists($excludes_fail_one{$variant})) {
            $excludes_fail_one{$variant} = 1;
        } else {
            # irgendwie ist die bedingung hier komisch. wenn die variante *nicht* im dict vorkommt,
            # soll sie hinzugefügt werden. falls doch, dann ins andere dict.
            $excludes_fail_two{$variant} = 1;
        }
    }
}

foreach(keys %excludes_fail_two) {
    print $perbatchfh "$_\n";
}
close $perbatchfh;

#
## Dump respective FDR thresholds and counts
#

open my $perbatch_t_fh, '>', $perbatch_thresholds_fn or die($!);

print $perbatch_t_fh "FDRthres\tFail_1plusbatches\tFail_2plusbatches\n";
for(my $f=0; $f<@fdr_thresholds; $f++) {
    %excludes_fail_one = ();
    %excludes_fail_two = ();
    for(my $b=0; $b < $batchcount; $b++) {
        for(my $i=0; $i < $perbatch_fdr[$b]->[$f]; $i++) {
            my $variant = $hwe_per_batch[$b]->[$i]->[0];
            if(not exists($excludes_fail_one{$variant})) {
                $excludes_fail_one{$variant} = 1;
            } else {
                # irgendwie ist die bedingung hier komisch. wenn die variante *nicht* im dict vorkommt,
                # soll sie hinzugefügt werden. falls doch, dann ins andere dict.
                $excludes_fail_two{$variant} = 1;
            }
        }
    }

    print $perbatch_t_fh sprintf("%e\t%d\t%d\n",
        $fdr_thresholds[$f],
        scalar(keys(%excludes_fail_one)),
        scalar(keys(%excludes_fail_two))
    );
}
close $perbatch_t_fh;

