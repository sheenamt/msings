#!/usr/bin/perl
#use strict;
use warnings;
use Test::More tests => 5;

use Stdev;

# Verify module can be included via "use" pragma
BEGIN { use_ok('Stdev') };

# Verify module can be included via "require" pragma
require_ok( 'Stdev' );

# Test true 
my @test_array1 = (1, 2, 3);
my $argumentTestCall1 = &Stdev::stdev(\@test_array1);
ok( $argumentTestCall1 == 1 );

# Test false
my $argumentTestCall2 = &Stdev::stdev(\@test_array1);
ok( $argumentTestCall2 != 3 );

# Test edge case
my @test_array2 = (1.1, 2.1, 31.1);
my $argumentTestCall3 = &Stdev::stdev(\@test_array2);
like( $argumentTestCall3 , '/17/' );

