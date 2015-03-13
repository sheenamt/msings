#mSINGS: Detecting Microsatellite Instability By Next-Generation Sequencing

##This project provides a program with a command-line interface for detecting microsatellite instability by next-generation sequencing.  This project is free for use by Academic users.  Please see enclosed license agreement for details and conditions of use.

##Citation:

Salipante SJ, Scroggins SM, Hampel HL, Turner EH, Pritchard CC.  2014. Microsatellite Instability Detection by Next Generation Sequencing.  Clin Chem. 2014 Jun 30.  [Epub ahead of print]

###Authors:  
   Steve Salipante and
   Sheena Scroggins
   
###Dependencies:  
   Perl 5x  
   Varscan (tested with 2.3.6)  
   Samtools (tested with 0.1.18)  

###Installation:

Clone git repo
  git clone git@bitbucket.org:uwlabmed/msings.git  

<OR> 

Download and extract zipped repo manually onto your system (go to "Downloads" section, at left)
   
Change directory ("cd") to mSINGS installation folder on your system before executing scripts.

###Input files:  
   1) bam file - sample of interest aligned against reference genome, provided in bam format  
   2) fasta file alignment was created with - in other words, your reference genome.  Must be indexed with bwa for use with samtools program (see below)
   3) MSI bed file (see example under "doc/mSINGS_TCGA.bed") - specifies the locations of the microsatellite tracts of interest.  NOTE:  must be sorted numerically and must not have a header line (see below). 
   4) MSI baseline file (see example under "doc/mSINGS_TCGA.baseline")  - describes the average and standard deviation of number expected peaks, as calculated from an MSI negative population (blood samples or MSI negative tumors).  User generates this file with create_baseline.py (see below) 
   5) MSI interval file (see example under "doc/mSINGS_TCGA.intervals")  - file for internal program use.  User makes this using MSI_formatter.pl  
   
Please note that both your reference genome and bed files MUST follow the convention that chromosomes are numbered numerically, without a "chr" prefix.  [i.e. - "10" instead of "chr10"].

###Execution:

Note: When specifying inputs, always provide absolute paths to files.  

1 Before you begin, ensure that your bed file is properly formatted.  Delete any header line, if present, then sort the file numerically:

    sort -k1,1n -k2,2n <unsorted_msi_bed.bed> -o <sorted_msi_bed.bed>

2 If necessary, create a bwa index for your reference genome

     bwa index <ref.fa>


3 Create the interval file, provide absolute paths in variables. 

    cd msings/bin

    perl -Ilib MSI_formatter.pl <output_dir> <msi_interval_file_name> <msi_bedfile>


4 Next, run the program for each MSI negative specimen to include in the baseline file creation

    bash call_msi.sh <bam_file> \
    <reference_genome_fasta> \
    <varscan.jar> \
    <sample name> \
    <output directory> \
    <msi_bed_file> \
    <msi_interval_file>
    

5 Create baseline file from MSI negative specimens

    python create_baseline.py <path to analysis files> <baseline file name>


###Note: Steps 1-5 only need to be done once per assay/target data set. Files may be saved and re-used for subsequent analyses.


6 Run the analysis script for each specimen of interest.  Specimens may be saved to the same output folder for convenience / sample batching.

    cd msings/bin  

    bash call_msi.sh <bam_file> \
    <reference_genome_fasta> \
    <varscan.jar> \
    <sample name> \
    <output directory> \
    <msi_bed_file> \
    <msi_interval_file>
    

7 After running MSI analysis for each specimen of interest, complete MSI analysis for a folder of output files.


    python count_msi.py <msi_baseline_file> <path to analysis files> <final output/analysis file name>


###Tests:

      cd msings/bin 

      prove -Ilib t

Should yield the following message:

      t/Average.t ........ ok

      t/Calc_MSI_dist.t .. ok

      t/Stdev.t .......... ok

      All tests successful.

      Files=3, Tests=18,  0 wallclock secs ( 0.03 usr  0.00 sys +  0.04 cusr  0.02 csys =  0.09 CPU)

      Result: PASS


https://bitbucket.org/uwlabmed/msings