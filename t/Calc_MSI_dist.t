#!/usr/bin/perl
#use strict;
use warnings;

use Test::More qw(no_plan);

use Calc_MSI_dist;

# Verify module can be included via "use" pragma
BEGIN { use_ok('Calc_MSI_dist') };

# Verify module can be included via "require" pragma
require_ok( 'Calc_MSI_dist' );

my %number_polymorphisms;
$number_polymorphisms{-1} = 1;
my %fraction_polymorphisms;
$fraction_polymorphisms{-1}= 0.0025974025974026;
my $totalmutantcount = 1;
my $variant = "INS-1-T:1:1:23:1:0:1";
my $uniquealleles = 1; 
my $totalcount = 280;

my($NEWnumber_polymorphisms, $NEWfraction_polymorphisms, $NEWtotalmutantcount, $NEWuniquealleles) = &Calc_MSI_dist::calc_MSI_dist(\%number_polymorphisms, \%fraction_polymorphisms, $totalmutantcount, $variant, $uniquealleles, $totalcount);
		
#overwrite old variables with the updated data and continue on
%number_polymorphisms=%{$NEWnumber_polymorphisms};
%fraction_polymorphisms=%{$NEWfraction_polymorphisms};
$totalmutantcount=$NEWtotalmutantcount;				
$uniquealleles=$NEWuniquealleles;

ok( $uniquealleles == 2 );
ok ( $totalmutantcount == 2);
ok ( $number_polymorphisms{1} == 1 );
ok ( $number_polymorphisms{-1} == 1 );
ok ( $fraction_polymorphisms{1} = 0.00357142857142857 );
ok ( $fraction_polymorphisms{-1} = 0.0025974025974026 );

