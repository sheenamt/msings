#!/usr/bin/perl
use strict;
use warnings;

package Average;

$Average::VERSION = '0.1';
#=================== http://edwards.sdsu.edu/research/index.php/kate/302-calculating-the-average-and-standard-deviation
sub average
{
        my($data) = @_;
        if (not @$data) {
                die("Sub average Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

