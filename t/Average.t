#!/usr/bin/perl
#use strict;
use warnings;
use Test::More tests => 5;

use Average;

# Verify module can be included via "use" pragma
BEGIN { use_ok('Average') };

# Verify module can be included via "require" pragma
require_ok( 'Average' );

# Test true 
my @test_array1 = (1, 2, 3);
my $argumentTestCall1 = &Average::average(\@test_array1);
ok( $argumentTestCall1 == 2 );

# Test false
my $argumentTestCall2 = &Average::average(\@test_array1);
ok( $argumentTestCall2 != 3 );

# Test edge case
my @test_array2 = (1.1, 2.1, 3.1);
my $argumentTestCall3 = &Average::average(\@test_array2);
ok( $argumentTestCall3 == 2.1 );

