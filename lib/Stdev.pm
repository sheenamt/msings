#!/usr/bin/perl
use strict;
use warnings;

use Average;

package Stdev;


$Std_dev::VERSION = '0.1';
#=================== http://edwards.sdsu.edu/research/index.php/kate/302-calculating-the-average-and-standard-deviation


sub stdev
{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &Average::average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}
