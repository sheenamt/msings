#!/usr/bin/perl
use strict;
use warnings;

package Calc_MSI_dist;

$Calc_MSI_dist::VERSION = '0.1';

######## http://stackoverflow.com/questions/10135642/perl-modifying-hash-reference-in-subroutine
######### http://stackoverflow.com/questions/12754414/perl-modifying-hash-in-subroutine-global-symbol-requires-explicit-package-name


sub calc_MSI_dist
{
# &calc_MSI_dist(\%number_polymorphisms, \%fraction_polymorphisms, $totalmutantcount, $variant, $uniquealleles, $totalcount);

my ($number_polymorphisms_hashref, $fraction_polymorphisms_hashref, $totalmutantcount, $variant, $uniquealleles, $totalcount) = @_;

		#declare internal variables
			my $adjustedcount = 0;
			my @splittemp;
			my $length;
			my @splitreads;

#copy the hash reference locally so that we can actually work with it as a hash.
my %number_polymorphisms=%{$number_polymorphisms_hashref};
my %fraction_polymorphisms=%{$fraction_polymorphisms_hashref};


					if ($variant =~ /^DEL/ or $variant =~ /^INS/)
						{
						#count it as a unique allele
						$uniquealleles ++;

						@splittemp = split("-", $variant);

						#if a deletion, negate the length, otherwise, use the stated length
						if ($variant =~ /^DEL/){$length = (-1 * $splittemp[1]);}
						else {$length = $splittemp[1];}

						#extract # variant reads
						@splitreads = split(":", $variant);

						#calculate fraction reads by dividing # reads by the wild-type at this site.  Need a provision in case there are no wild type reads (unlikely!)
						if ($totalcount ne 0)
						{
						$adjustedcount = ($splitreads[1] / $totalcount);
						}
						else {$adjustedcount = 0;}

						#extract mutant count

						$totalmutantcount = $totalmutantcount + $splitreads[1];
#	print "$startposition	$variant	$splitreads[1]	$totalmutantcount\n";
						#check to see if this value is already in our hash
						#if this indel has been seen before, increase the value of the hash by the fraction of reads
						if (exists $fraction_polymorphisms{$length}){$fraction_polymorphisms{$length} = $fraction_polymorphisms{$length} + $adjustedcount;}
						#if we have not seen this indel before, add it to the hash with a value of one
						else {$fraction_polymorphisms{$length} = $adjustedcount;}

						if (exists $number_polymorphisms{$length}){$number_polymorphisms{$length} = $number_polymorphisms{$length} + $splitreads[1];}
						else {$number_polymorphisms{$length} = $splitreads[1];}

						}

							
						
				return (\%number_polymorphisms, \%fraction_polymorphisms, $totalmutantcount, $uniquealleles);
						

					
}
