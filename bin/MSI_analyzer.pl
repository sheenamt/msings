#!/usr/bin/perl

use Stdev;
use Calc_MSI_dist;

#OPERATIONS:
# reads in specified bed file, counts #, percentage, and distribution of indels (not base substitutions) within each range of the bed
#incorporates cutoff for minimum peak height - set to 2.5%.  corrects some glitches with off-target counting of reads.

#accept arguments from command line
my $savepath = $ARGV[0];
my $prefix = $ARGV[1];
my $msi_bed_file = $ARGV[2]; #include path to bed file to use

#my $savepath = "./";
#my $prefix = "test";#"LMG-240";#
# my $msi_bed_file = "hg19_comprehensive_50bp_V3.bed";


my $peakfractioncutoff = 0.05;

my $datafile = $prefix.".msi_output";


#prepare file for output
my $segmentedoutputfile = $prefix.".msi.txt";
	open( output, ">$savepath/$segmentedoutputfile");
	print output "Position	Name	Avg_read_depth	Number_Peaks	Standard_Deviation	Peak_distribution-->\n";
	close( output);

#arrays to save average values of stdev and and peak numbers
my @avestdevarray;
my @avepeaknumberarray;

#open the list of msi sites, consider each in turn
my $bedfileFH = "CAP";
my $datafileFH = "MSI";
open ($bedfileFH, "$msi_bed_file");
open ($datafileFH, "$savepath/$datafile");

#check to be sure that the ouput we are readin in is not empty
my $linecount = 0;

while(<MSI>)
{
    $linecount ++;
    if ($linecount > 2) {last;}
}

if ($linecount < 2){die "input file $datafile is empty.  Aborting!\n" ;}


sub read_file_line {
  my $fh = shift;

  if ($fh and my $line = <$fh>) {
    chomp $line;
    return $line;
  }
  return;
}

	#read initial lines from both files
	
my $line = &read_file_line($bedfileFH);

#read twice to skip header!!
my $sline = &read_file_line($datafileFH);
$sline = &read_file_line($datafileFH);

#for each line in the BED file do the following.  initialize here.
	#prepare a hash to keep track of each size of indel and the fraction of reads per each
	my %fraction_polymorphisms;
	#prepare a hash to keep track of each size of indel and the number of reads per each
	my %number_polymorphisms;


	#prepare a counter to keep track of the number of unique alleles
	my $uniquealleles = 0;

	#counter to keep track of cumulative read depth
	my $readdepth = 0;

	#and number of bases considered
	my $totalbases = 0;
	#count the number of mutants reads we encounter.  WT will by total reads - mutant!
	my $totalmutantcount = 0;

while ($bedfileFH ) #go through whole bed file
	{

	#split out lines from file to an array;
	my @site = split ("\t", $sline);
	my @msi = split ("\t", $line);

	my $chromosome = $msi[0];
	my $startposition = $msi[1];
	my $stopposition = $msi[2];

	my $currentposition;

			$startposition =~ s/\s//g;

			$stopposition =~ s/\s//g;

		if (eof($bedfileFH))
			{
					open( output, ">>$savepath/$segmentedoutputfile");
					print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	0	1	0	0:0:0";

					print output "\n";

					close( output);			
			}
			
		last if (eof($bedfileFH));
			
			
			#first tackle reasons to advance the BED file
			#1)if the bed file chromosome is less than the chromosome in the datafile, advance the line in the bed file	
			if ($chromosome < $site[0])
				{ #print "contition 1 \n";				
				
					#Calculate summary statistics and print to file

					#average read depth is totalcount divided by # bases considered
					my $averagedepth = 0;
					if ($totalbases ne 0){$averagedepth = int($readdepth / $totalbases);}

					#calculate fraction and number of wild-type alleles
				print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	\n";
					my $wildtypefraction = 0;
					if ($averagedepth ne 0 and $totalmutantcount < $averagedepth){$wildtypefraction = ($averagedepth - $totalmutantcount) / $averagedepth;}
						elsif ($averagedepth ne 0 and $totalmutantcount >= $averagedepth){$wildtypefraction = 0;}
					my $wildtypenumber = int($averagedepth - $totalmutantcount);


					#get the max and min key from the hash
					my $minhash = 0;
					my $maxhash = 0;
					my $numberpeaks = 1;

					while(($key,$value) = each(%fraction_polymorphisms))
						{
						if ($value >= $peakfractioncutoff ) {$numberpeaks ++;}
						if ($key < $minhash){$minhash = $key;}
						if ($key > $maxhash){$maxhash = $key;}
						}

					#make an array with each value in this range
					my @variants = ($minhash .. $maxhash);


					#calculate standard deviation of alleles
					my @stdevarray;
					while(($key,$value) = each(%number_polymorphisms))
						{
						my $count = 0;
						while ($count < $value)
							{
							push(@stdevarray, $key);
							$count ++;
							}
						}
					#add in wt
					my $count2=0;
					while ($count2 < $wildtypenumber)
						{
						push(@stdevarray,"0");
						$count2++;
						}

					if (scalar @stdevarray >= 2)
					{
					$std = &Stdev::stdev(\@stdevarray);
					}
					else { $std = 0;}

					open( output, ">>$savepath/$segmentedoutputfile");
					print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	$std";

					push(@avestdevarray, $std);
					push(@avepeaknumberarray, $numberpeaks);

					#now print out the keys and values for each item in the hash
					foreach my $polymorphism (@variants)
						{
						#special provision for zero, which is wild type polymorphism
						if ($polymorphism eq "0")
							{
							print output "	0:","$wildtypefraction",":","$wildtypenumber";
							}
						elsif (exists $fraction_polymorphisms{$polymorphism})
							{
							print output "	$polymorphism",":","$fraction_polymorphisms{$polymorphism}",":","$number_polymorphisms{$polymorphism}";
							}
						else
							{
							print output "	$polymorphism",":","0";
							}

						}

					print output "\n";

					close( output);
	
				$line = &read_file_line($bedfileFH);
	
				#prepare a hash to keep track of each size of indel and the fraction of reads per each
				%fraction_polymorphisms = ();
				#prepare a hash to keep track of each size of indel and the number of reads per each
				%number_polymorphisms = ();

				#prepare a counter to keep track of the number of unique alleles
				$uniquealleles = 0;

				#counter to keep track of cumulative read depth
				$readdepth = 0;

				#and number of bases considered
				$totalbases = 0;
				#count the number of mutants reads we encounter.  WT will by total reads - mutant!
				$totalmutantcount = 0;
	
				}
			#2)chromosome matches, but position in bed file is below that in the data file, advance the line in the bed file	
			elsif ($chromosome == $site[0] and $startposition < $site[1] and $stopposition < $site[1])
				{ #print "contition 2 \n";
				
					#Calculate summary statistics and print to file

					#average read depth is totalcount divided by # bases considered
				        my $averagedepth = 0;
					if ($totalbases ne 0){$averagedepth = int($readdepth / $totalbases);}

					#calculate fraction and number of wild-type alleles
				print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	\n";
					my $wildtypefraction = 0;
					if ($averagedepth ne 0 and $totalmutantcount < $averagedepth){$wildtypefraction = ($averagedepth - $totalmutantcount) / $averagedepth;}
						elsif ($averagedepth ne 0 and $totalmutantcount >= $averagedepth){$wildtypefraction = 0;}
					my $wildtypenumber = int($averagedepth - $totalmutantcount);


					#get the max and min key from the hash
					my $minhash = 0;
					my $maxhash = 0;
					my $numberpeaks = 1;

					while(($key,$value) = each(%fraction_polymorphisms))
						{
						if ($value >= $peakfractioncutoff ) {$numberpeaks ++;}
						if ($key < $minhash){$minhash = $key;}
						if ($key > $maxhash){$maxhash = $key;}
						}

					#make an array with each value in this range
					my @variants = ($minhash .. $maxhash);


					#calculate standard deviation of alleles
					my @stdevarray;
					while(($key,$value) = each(%number_polymorphisms))
						{
						my $count = 0;
						while ($count < $value)
							{
							push(@stdevarray, $key);
							$count ++;
							}
						}
					#add in wt
					my $count2=0;
					while ($count2 < $wildtypenumber)
						{
						push(@stdevarray,"0");
						$count2++;
						}

					if (scalar @stdevarray >= 2)
					{
					$std = &Stdev::stdev(\@stdevarray);
					}
					else { $std = 0;}

					open( output, ">>$savepath/$segmentedoutputfile");
					print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	$std";

					push(@avestdevarray, $std);
					push(@avepeaknumberarray, $numberpeaks);

					#now print out the keys and values for each item in the hash
					foreach my $polymorphism (@variants)
						{
						#special provision for zero, which is wild type polymorphism
						if ($polymorphism eq "0")
							{
							print output "	0:","$wildtypefraction",":","$wildtypenumber";
							}
						elsif (exists $fraction_polymorphisms{$polymorphism})
							{
							print output "	$polymorphism",":","$fraction_polymorphisms{$polymorphism}",":","$number_polymorphisms{$polymorphism}";
							}
						else
							{
							print output "	$polymorphism",":","0";
							}

						}

					print output "\n";

					close( output);
	
				$line = &read_file_line($bedfileFH);
	
				#prepare a hash to keep track of each size of indel and the fraction of reads per each
				%fraction_polymorphisms = ();
				#prepare a hash to keep track of each size of indel and the number of reads per each
				%number_polymorphisms= ();


				#prepare a counter to keep track of the number of unique alleles
				$uniquealleles = 0;

				#counter to keep track of cumulative read depth
				$readdepth = 0;

				#and number of bases considered
				$totalbases = 0;
				#count the number of mutants reads we encounter.  WT will by total reads - mutant!
				$totalmutantcount = 0;
	
				}				
			#3)chromosome matches, but position in bed file is higher than that in the data file, advance the line in the bed file	
			elsif ($chromosome == $site[0] and $startposition > $site[1] and $stopposition > $site[1])
				{ #print "contition 3 \n";
				
					#Calculate summary statistics and print to file

					#average read depth is totalcount divided by # bases considered
					my $averagedepth = 0;
					if ($totalbases ne 0){$averagedepth = int($readdepth / $totalbases);}

					#calculate fraction and number of wild-type alleles
				print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	\n";
					my $wildtypefraction = 0;

					if ($averagedepth ne 0 and $totalmutantcount < $averagedepth){$wildtypefraction = ($averagedepth - $totalmutantcount) / $averagedepth;}
						elsif ($averagedepth ne 0 and $totalmutantcount >= $averagedepth){$wildtypefraction = 0;}
					my $wildtypenumber = int($averagedepth - $totalmutantcount);



					#get the max and min key from the hash
					my $minhash = 0;
					my $maxhash = 0;
					my $numberpeaks = 1;

					while(($key,$value) = each(%fraction_polymorphisms))
						{
						if ($value >= $peakfractioncutoff ) {$numberpeaks ++;}
						if ($key < $minhash){$minhash = $key;}
						if ($key > $maxhash){$maxhash = $key;}
						}

					#make an array with each value in this range
					my @variants = ($minhash .. $maxhash);


					#calculate standard deviation of alleles
					my @stdevarray;
					while(($key,$value) = each(%number_polymorphisms))
						{
						my $count = 0;
						while ($count < $value)
							{
							push(@stdevarray, $key);
							$count ++;
							}
						}
					#add in wt
					my $count2=0;
					while ($count2 < $wildtypenumber)
						{
						push(@stdevarray,"0");
						$count2++;
						}

					if (scalar @stdevarray >= 2)
					{
					$std = &Stdev::stdev(\@stdevarray);
					}
					else { $std = 0;}

					open( output, ">>$savepath/$segmentedoutputfile");
					print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	$std";

					push(@avestdevarray, $std);
					push(@avepeaknumberarray, $numberpeaks);

					#now print out the keys and values for each item in the hash
					foreach my $polymorphism (@variants)
						{
						#special provision for zero, which is wild type polymorphism
						if ($polymorphism eq "0")
							{
							print output "	0:","$wildtypefraction",":","$wildtypenumber";
							}
						elsif (exists $fraction_polymorphisms{$polymorphism})
							{
							print output "	$polymorphism",":","$fraction_polymorphisms{$polymorphism}",":","$number_polymorphisms{$polymorphism}";
							}
						else
							{
							print output "	$polymorphism",":","0";
							}

						}

					print output "\n";

					close( output);

	
				$line = &read_file_line($bedfileFH);
	
				#prepare a hash to keep track of each size of indel and the fraction of reads per each
				%fraction_polymorphisms = ();
				#prepare a hash to keep track of each size of indel and the number of reads per each
				%number_polymorphisms = ();


				#prepare a counter to keep track of the number of unique alleles
				$uniquealleles = 0;

				#counter to keep track of cumulative read depth
				$readdepth = 0;

				#and number of bases considered
				$totalbases = 0;
				#count the number of mutants reads we encounter.  WT will by total reads - mutant!
				$totalmutantcount = 0;
	
				}
	




			#now, reasons to advance DATA FILE
	
			#if the data file chromosome is less than the chromosome in the bed file, advance the line in the data file	
			elsif ($chromosome > $site[0])
				{
				$sline = &read_file_line($datafileFH);
			
				}
			#chromosome matches, but position in data file is below that in the bed file, advance the line in the data file	
			elsif ($chromosome == $site[0] and $site[1] < $startposition )
				{
				$sline = &read_file_line($datafileFH);
				}		

				
				
				
			#proceed if this site is contained within an msi site
			elsif ( $site[0] eq $chromosome and $site[1] >= $startposition and $site[1] <= $stopposition)
				{


				my $totalcount;
				
				#extract the total count of mapped reads
				$totalcount = $site[4];
				$readdepth = $readdepth + $totalcount;

				#add 1 to the total bases counted
				$totalbases ++;

				#consider alternative variants, referring to insertions and deletions separately
				foreach my $variant (@site)
					{

				#pass to subroutine, return new values for the hashes and appropriate variables.
				my($NEWnumber_polymorphisms, $NEWfraction_polymorphisms, $NEWtotalmutantcount, $NEWuniquealleles) = &Calc_MSI_dist::calc_MSI_dist(\%number_polymorphisms, \%fraction_polymorphisms, $totalmutantcount, $variant, $uniquealleles, $totalcount);

			
				#overwrite old variables with the updated data and continue on
				%number_polymorphisms=%{$NEWnumber_polymorphisms};
				%fraction_polymorphisms=%{$NEWfraction_polymorphisms};
				$totalmutantcount=$NEWtotalmutantcount;				
				$uniquealleles=$NEWuniquealleles;
				
			

					}
				#lastly, advance data file !
				$sline = &read_file_line($datafileFH);

				}

				# final special provision : if we have hit the end of the data file  need to print the very last entry, if non-empty:
				if ($sline eq "") #######################
				{
				unless (!keys %fraction_polymorphisms)
				{
					#Calculate summary statistics and print to file
					#average read depth is totalcount divided by # bases considered
					my $averagedepth = 0;
					if ($totalbases ne 0){$averagedepth = int($readdepth / $totalbases);}
					
					#calculate fraction and number of wild-type alleles
				print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	\n";
					my $wildtypefraction = 0;
					if ($averagedepth ne 0 and $totalmutantcount < $averagedepth){$wildtypefraction = ($averagedepth - $totalmutantcount) / $averagedepth;}
						elsif ($averagedepth ne 0 and $totalmutantcount >= $averagedepth){$wildtypefraction = 0;}
					my $wildtypenumber = int($averagedepth - $totalmutantcount);

				#print " avedepth: $averagedepth	totalmutant: $totalmutantcount	wiltypefrac: $wildtypefraction wildtypenumber:	$wildtypenumber\n\n";

					#get the max and min key from the hash
					my $minhash = 0;
					my $maxhash = 0;
					my $numberpeaks = 1;

					while(($key,$value) = each(%fraction_polymorphisms))
						{
						if ($value >= $peakfractioncutoff ) {$numberpeaks ++;}
						if ($key < $minhash){$minhash = $key;}
						if ($key > $maxhash){$maxhash = $key;}
						}

					#make an array with each value in this range
					my @variants = ($minhash .. $maxhash);


					#calculate standard deviation of alleles
					my @stdevarray;
					while(($key,$value) = each(%number_polymorphisms))
						{
						my $count = 0;
						while ($count < $value)
							{
							push(@stdevarray, $key);
							$count ++;
							}
						}
					#add in wt
					my $count2=0;
					while ($count2 < $wildtypenumber)
						{
						push(@stdevarray,"0");
						$count2++;
						}

					if (scalar @stdevarray >= 2)
					{
					$std = &Stdev::stdev(\@stdevarray);
					}
					else { $std = 0;}

					open( output, ">>$savepath/$segmentedoutputfile");
					print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	$averagedepth	$numberpeaks	$std";

					push(@avestdevarray, $std);
					push(@avepeaknumberarray, $numberpeaks);

					#now print out the keys and values for each item in the hash
					foreach my $polymorphism (@variants)
						{
						#special provision for zero, which is wild type polymorphism
						if ($polymorphism eq "0")
							{
							print output "	0:","$wildtypefraction",":","$wildtypenumber";
							}
						elsif (exists $fraction_polymorphisms{$polymorphism})
							{
							print output "	$polymorphism",":","$fraction_polymorphisms{$polymorphism}",":","$number_polymorphisms{$polymorphism}";
							}
						else
							{
							print output "	$polymorphism",":","0";
							}

						}

					print output "\n";

					close( output);	

					#advance the bed line.
					$line = &read_file_line($bedfileFH);
					#prepare a hash to keep track of each size of indel and the fraction of reads per each
					%fraction_polymorphisms = ();
					#prepare a hash to keep track of each size of indel and the number of reads per each
					%number_polymorphisms = ();

					#prepare a counter to keep track of the number of unique alleles
					$uniquealleles = 0;

						#counter to keep track of cumulative read depth
					$readdepth = 0;
	
					#and number of bases considered
					$totalbases = 0;
					#count the number of mutants reads we encounter.  WT will by total reads - mutant!
					$totalmutantcount = 0;

					}
				else
					{
				#if we have reached the end of the data file, still need to print the bed file line with zero depth reported.  ###########

					open( output, ">>$savepath/$segmentedoutputfile");
					print output "$chromosome",":","$startposition","-","$stopposition	$msi[3]	0	1	0	0:0:0";

					print output "\n";

					close( output);
	
					#advance the bed line.
					$line = &read_file_line($bedfileFH);
					#prepare a hash to keep track of each size of indel and the fraction of reads per each
					%fraction_polymorphisms = ();
					#prepare a hash to keep track of each size of indel and the number of reads per each
					%number_polymorphisms = ();

					#prepare a counter to keep track of the number of unique alleles
					$uniquealleles = 0;

						#counter to keep track of cumulative read depth
					$readdepth = 0;
	
					#and number of bases considered
					$totalbases = 0;
					#count the number of mutants reads we encounter.  WT will by total reads - mutant!
					$totalmutantcount = 0;

									
					}
				}



	}

	
exit;



