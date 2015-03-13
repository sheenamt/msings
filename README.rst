mSINGS: Detecting Microsatellite Instability By Next-Generation Sequencing
============================================================================
This project provides a program with a command-line interface for detecting microsatellite instability by next-generation sequencing.  This project is free for use by Academic users.  Please see enclosed license agreement for details and conditions of use.

Citation:

Salipante SJ, Scroggins SM, Hampel HL, Turner EH, Pritchard CC.  2014. Microsatellite Instability Detection by Next Generation Sequencing.  Clin Chem. 2014 Sep;60(9):1192-9. doi: 10.1373/clinchem.2014.223677. Epub 2014 Jun 30.

Authors:

* Steve Salipante
* Sheena Scroggins

Dependencies:
-------------

* Tested on Ubuntu 12.04
* Python 2.7.x
* git
* Installed in virtualenv:

 * SAMTOOLS (VERSION=0.1.18)
 * VARSCAN (VERSION=v2.3.7)

Installation:
-------------
Clone repo using git
     $ git clone git@bitbucket.org:uwlabmed/msings.git

Or download and extract zipped repo manually onto your system (go to "Downloads" section, at left)

Change directory ("cd") to mSINGS installation folder on your system before executing scripts.
     $ cd path/to/msings

Creating a virtualenv for running pipeline:
--------------------------------------------
This pipeline is designed to be run inside its own virtualenv, specified in settings.conf. To create a new virtualenv, first navigate to the top level directory of the repo and run the install script, specifying a directory for the virtualenv if pwd-env isn't sufficient:
     $ cd path/to/msings

     $ bash dev/bootstrap.sh /different/virtual/path

To use the virtualenv, including the programs (varscan, samtools, msi), run:
     $ source msings-env/bin/activate

Development (developers only):
-----------
The devel.sh script builds a local virtualenv and downloads test data (if run from UW):
     $ git clone git@bitbucket.org:uwlabmed/msings.git

     $ bash dev/devel.sh 

Required Input files:
------------
1. bam file - sample of interest aligned against reference genome, provided in bam format. Index required. 

2. ref_fasta : fasta file alignment was created with - in other words, your reference genome.  Must be indexed with bwa for use with samtools program (see below).  Please note that both your reference genome and bed files MUST follow the convention that chromosomes are numbered numerically or with "X" or "Y".  Other conventions (such as those bearing a "chr" prefix) are not supported.

3. msi_bed : MSI bed file (see example under "doc/mSINGS_TCGA.bed") - specifies the locations of the microsatellite tracts of interest.  NOTE:  must be sorted numerically and must not have a header line (see below), and must follow identical naming conventions as the reference genome.

4. msi_baseline : MSI baseline file (see example under "doc/mSINGS_TCGA.baseline")  - describes the average and standard deviation of the number of expected signal peaks at each locus, as calculated from an MSI negative population (blood samples or MSI negative tumors).  User generates this file with msi create_baseline script (see below).

5. msi_intervals : MSI interval file (see example under "doc/mSINGS_TCGA.intervals")  - file for internal program use.  User makes this using msi formatter script.

6. settings.conf : file specifying absolute paths to above files, thresholds for MSI status calling, and control name.

7. data.conf : file specifying absolute paths to bam files to be processed. 

Please note that both your reference genome and bed files MUST follow the convention that chromosomes are numbered numerically or with "X" or "Y" - other names are not supported.


Output files:
------------
For each sample run, the following will be produced:
 * SAMPLE.MSI.txt = very detailed information about instability and distribution of alleles of differing length.  Raw data that is used to generate final MSI calls.
 * SAMPLE.MSI_output.txt = output of Varscan readcounts function, an intermediate file.
 * SAMPLE.MSI_Analysis.txt = Binary matrix of interpreted instability (1) or stability (0) at each locus. Loci with insufficient coverage for instability calling are left blank. Summary statistics and interpretation of results are provided.

For the entire run, a "top level" output represented as a binary matrix of interpreted instability (1) or stability (0) at each locus is provided if the count_msi.py function is run. Loci with insuffient coverage for instability calling are left blank. Summary statistics and interpretation of results are provided.

Execution for Exome / TCGA data sets:
----------
This protocol will run the pipeline using the baseline file and microsatellite loci identified for TCGA exome data.
Files specific for analysis of TCGA exome data are provided in the doc/ directory of this package. 
 * msi_baseline 
 * msi_bed 
 * msi_intervals 

1. Edit the settings.conf to point to the absolute path of the input bams, the name of the control sample, and the path to the reference fasta used to align this bam:

    datadir = /path/to/inputbamfiles  
    baseline_control = control_name
    ref_fasta = /path/to/ref.fasta

2. To test that everything is installed and all inputs are specified correctly, the -n flag can be used: 

  $ scons -n
  
  scons: Building targets

  <......>

  scons: done building targets.

3. Run the analysis script for the batch of samples. Output will be in the output directory specified in the settings.conf file, 'output' by default

 $ scons 

Execution for custom data sets:
-------------------
Files specific for analysis of TCGA exome data are provided in the doc/ directory of this package. To run mSINGS analysis use custom assays or custom targets, users are required to provide 3 custom files:
 * msi_baseline 
 * msi_bed 
 * msi_intervals 

The folllwing instructions will allow users to set up analysis for their custom targets, to generate a custom baseline for those targets, and to run supsequent analysis.

1. Before you begin creating custom files, acivate the virtualenv to make use of installed programs:
  
 $ source msings-env/bin/activate

2. Ensure that your bed file is properly formatted.  Delete any header line, if present, then sort the file numerically by chromosome and then by base position:

 $ sort -V -k1,1 -k2,2n /path/to/CUSTOM_UNSORTED_MSI_BED -o /path/to/CUSTOM_MSI_BED

3. If necessary, bwa format and create a bwa index for your reference genome:

 $ bwa index -a bwtsw ref_fasta

4. Create the interval file, providing absolute paths for variables:

 $ msi formatter /path/to/CUSTOM_MSI_BED -o /path/to/CUSTOM_MSI_INTERVALS

5. Now that we have CUSTOM_MSI_BED and CUSTOM_MSI_INTERVALS, update the settings.conf to reflect these:

  msi_bed = /path/to/custom_MSI_BED

  msi_intervals = /path/to/CUSTOM_MSI_INTERVALS

6. Run analysis to generate a baseline file from MSI negative specimens.  Edit the settings.conf file to point to the absolute path of the MSI negative specimen BAMS and provide a name for the desired output folder: 

  #Absolute path to input bam files for MSI negative specimens to be included in the baseline

  datadir = /path/to/inputbamfiles 

  #output folder name

  output = /path/to/my_output


7. Now test the setup for the creation of the msi-calls files, and run the program for each MSI negative specimen to include in the baseline file creation:

 $ scons -n msi-calls

 If that produces the expected "scons: done building targets" message, proceed: 

 $ scons msi-calls

8. Use these raw data to produce the the MSI_BASELINE file from MSI negative specimens

 $ msi create_baseline /path/to/my_output -o /path/to/CUSTOM_MSI_BASELINE

NOTE: This process only need to be done once per assay/target data set. Files may be saved and re-used for subsequent analyses. 

9. Now we update the settings.conf to point to all the new custom files:

 msi_bed = /path/to/CUSTOM_MSI_BED

 msi_intervals = /path/to/CUSTOM_MSI_INTERVALS

 msi_baseline = /path/to/CUSTOM_MSI_BASELINE
 
Once the settings.conf file is updated for the new custom files, execution is the same as for Exome / TCGA data sets (above).  

10. To test that everything is installed and all inputs are specified correctly, the -n flag can be used: 

  $ scons -n
  
  scons: Building targets

  <......>

  scons: done building targets.

11. Run the analysis script for the batch of samples. Output will be in the output directory specified in the settings.conf file, 'output' by default

 $ scons 
 
Tests:
^^^^^^

   cd msings
   source msings-env/bin/active
   ./testall

Should yield the following message:
 Ran 11 tests in 0.068s

 OK

https://bitbucket.org/uwlabmed/msings