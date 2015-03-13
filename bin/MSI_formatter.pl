#!/usr/bin/perl

#OPERATIONS:
#formats MSI bed file for input into varscan's readcount function.
# NOTE - bed file only requires first 3 lines. Include the rest just to keep track.  regions include tract of interest and 3-5 bases of flanking

#only care about the positions, not the variant, so just consider everything as an insertion
#about format
#Sure, thanks for asking. I should tell you that the "readcounts" function is almost deprecated, as the methods for consensus calling (e.g. mpileup2cns) are more comprehensive and robust. That said, if you give the 'readcounts' command a list of variant positions (--variant-file), the file should be in 4-column tab delimited format: 
#chromosome 
#position 
#reference base 
#variant base

#accept arguments from command line
my $savepath = $ARGV[0];
my $prefix = $ARGV[1];
my $msi_file = $ARGV[2];

#my $savepath = "./";
#my $prefix = "test";
#my $msi_file = "OncoPlex_MSI_v3_03052013.txt";

#prepare file for output
my $segmentedoutputfile = $prefix.".msi_intervals";
	open( output, ">$savepath/$segmentedoutputfile");
	close( output);	

#open the list of clinical variants, exclude header line
open (CAP, "$msi_file");
while (<CAP>) 
	{
	my $line = $_;
	chomp $line;
	
	my @site;
	my $chromosome;
	my $startposition;
	my $stopposition;
	
	my $currentposition;

		#split to components
		@site = split ("\t", $line);
		
		
		#do nOT! remove "chr" prefix
		$chromosome = $site[0];
		$chromosome =~ s/chr//g;
				
		$startposition = $site[1];
			$startposition =~ s/\s//g;
		$stopposition = $site[2];
			$stopposition =~ s/\s//g;
		
	#make an entry for each position in the range
			
		$currentposition = $startposition;
		
			while ($currentposition <= $stopposition)
			{
			
			#now print to output in tab-delimited txt.
			open( output, ">>$savepath/$segmentedoutputfile");
			print output "$chromosome	$currentposition	-	T\n";
			close( output);	
			
			$currentposition ++;
			
			}

		
	}

exit;
