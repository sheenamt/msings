==========================================================================
mSINGS: Detecting Microsatellite Instability By Next-Generation Sequencing
==========================================================================

This project provides a program with a command-line interface for detecting microsatellite instability by next-generation sequencing.  This project is free for use by Academic users.  Please see enclosed license agreement for details and conditions of use.


Citation:
^^^^^^^^^
Salipante SJ, Scroggins SM, Hampel HL, Turner EH, Pritchard CC.  2014. Microsatellite Instability Detection by Next Generation Sequencing.  Clin Chem. 2014 Sep;60(9):1192-9. doi: 10.1373/clinchem.2014.223677. Epub 2014 Jun 30.

Authors:
^^^^^^^^
* Steve Salipante
* Sheena Todhunter


Dependencies:
^^^^^^^^^^^^^
* Tested on Ubuntu 12.04, 14.04
* Python 2.7.x
* git
* virtualenv

Installation:
-------------
Please ensure the dependencies are installed to the system. 
This pipeline is designed to be run inside its own virtual environment (virtualenv) which contains all required dependencies. The following commands will allow you to establish this virtualenv.

To install the software to a virtual environment in your current directory, run the following commands:

>>> git clone https://bitbucket.org/uwlabmed/msings.git
>>> cd msings
>>> bash dev/bootstrap.sh

To use the virtualenv, including the programs (varscan, samtools, msi), run:

>>> source msings-env/bin/activate

To install the software (including creating a virtualenv) to a path of your choosing, run the following commands:

>>> git clone https://bitbucket.org/uwlabmed/msings.git
>>> cd msings
>>> bash dev/bootstrap.sh /desired/virtual/path

NOTE:  Downloading the zipped file through this website will cause installation problems and is not supported. Please CLONE the repo. 

Development (UW developers only):
------------------------------
The devel.sh script builds a local virtualenv and downloads test data (if run from UW):

>>> git clone https://bitbucket.org/uwlabmed/msings.git
>>> cd msings 
>>> bash dev/devel.sh 

If run outside of the UW network, it will not download test data. 

Required Input files:
----------------------
1. bam file : sample of interest aligned against reference genome, provided in bam format. Index required. 

2. ref_fasta : fasta file alignment was created with - in other words, your reference genome.  Must be indexed with bwa for use with samtools program (see below).  Please note that both your reference genome and bed files MUST follow the convention that chromosomes are numbered numerically or with "X" or "Y".  Other conventions (such as those bearing a "chr" prefix) are not supported.

3. msi_bed : MSI bed file (see example under "doc/mSINGS_TCGA.bed") - specifies the locations of the microsatellite tracts of interest.  NOTE:  must be sorted numerically and must not have a header line (see below), and must follow identical chromosome naming conventions as the reference genome.

4. msi_baseline : MSI baseline file (see example under "doc/mSINGS_TCGA.baseline")  - describes the average and standard deviation of the number of expected signal peaks at each locus, as calculated from an MSI negative population (blood samples or MSI negative tumors).  User generates this file with msi create_baseline script (see below).  **IMPORTANT NOTE**:  Baseline statistics vary markedly from assay-to-assay and lab-to-lab.  It is **CRITICAL** that you prepare a baseline file that is specific for your analytic process, and for which data have been generated using the same protocols. 

5. msi_intervals : MSI interval file (see example under "doc/mSINGS_TCGA.intervals")  - file for internal program use.  User makes this using create_intervals script.

Please note that both your reference genome and bed files MUST follow the same naming convention of chromosomes (ie "chr1" or just "1", without the "chr" prefix).  Chromosomes are numbered numerically, or as "X" or "Y" - other chromsome names are not supported.

Output files:
-------------
For each sample run, the following will be produced:
 * SAMPLE.msi.txt = very detailed information about instability and distribution of alleles of differing length.  Raw data that is used to generate final MSI calls.
 * SAMPLE.msi_output.txt = output of Varscan readcounts function, an intermediate file.
 * SAMPLE.MSI_Analysis.txt = Binary matrix of interpreted instability (1) or stability (0) at each locus. Loci with insufficient coverage for instability calling are left blank. Summary statistics and interpretation of results are provided.

For the entire run, a "top level" output represented as a binary matrix of interpreted instability (1) or stability (0) at each locus is provided if the count_msi.py function is run. Loci with insufficient coverage for instability calling are left blank. Summary statistics and interpretation of results are provided.

Execution for Exome / TCGA data sets: 
-------------------------------------
This protocol will run the pipeline using the baseline file and microsatellite loci identified for TCGA exome data. Files specific for analysis of TCGA exome data are provided in the doc/ directory of this package. 
 * msi_baseline 
 * msi_bed 
 * msi_intervals 

1. If you installed the virtualenv to a different location that the default scripts, you MUST edit the bash scripts to point to your virtual environment and your VarScan jar file

>>> source /path/to/your/msings-virtual-environment/bin/activate

2. OPTIONAL - Edit the run_msings.sh to point to the absolute path of the reference genome used for alignment. If you choose to not edit the script, you will be required to specify the reference genome file at the command line in order to execute the script

>>> REF_GENOME=/path/to/REF_GENOME;

3. OPTIONAL - Edit the run_msings.sh to change the MSI default analytic parameters:
 
>>> multiplier=2.0 
    "multiplier" is the number of standard deviations from the baseline that is required to call instability
   
>>> msi_min_threshold=0.2
    "msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative   

>>> msi_max_threshold=0.2
    "msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive

* If the fraction of unstable sites falls between the thresholds, the specimen is considered indeterminate.  (By default, no indeterminate calls are permitted) 

4. Create a file of the list of BAMS, with each line being the absolute path to one sample

>>> /path/to/sampleA.bam
>>> /path/to/sampleB.bam
>>> /path/to/sampleC.bam
   
5. Run the analysis script for the batch of samples. Output will be in subfolders of the same directories containing the BAM data, subfolders are named after the samples themselves

Default execution:
>>>  scripts/run_msings.sh PATH/TO/BAM_LIST PATH/TO/INTERVALS_FILE PATH/TO/BEDFILE PATH/TO/REF_GENOME PATH/TO/MSI_BASELINE

If you already edited the run_msings.sh script to point to the reference files (either yours or the TCGA files in the doc/ folder), script may be run as follows:
>>>  scripts/run_msings.sh PATH/TO/BAM_LIST

Execution for custom data sets:
-------------------
Files specific for analysis of TCGA exome data are provided in the doc/ directory of this package. To run mSINGS analysis use custom assays or custom targets, users are required to provide 3 custom files:
 * msi_baseline 
 * msi_bed 
 * msi_intervals 

NOTE: loci PRESENT in the bed file that are ABSENT in the baseline file (created in step 8 below) will not be scored!

The following instructions will allow users to set up analysis for their custom targets, to generate a custom baseline for those targets, and to run subsequent analysis.  Recommendations for design of custom assays and custom targets are provided in the Recommendations_for_custom_assays.txt file packaged with the repository.

1. If you installed the virtualenv to a different location than the default scripts, you MUST edit the following bash scripts to point to your specific virtual environment and VarScan jar file

scripts/create_intervals.sh:
>>> source /path/to/your/msings-virtual-environment/bin/activate
>>> VARSCAN=/path/to/your/msings-virtual-environment/bin/VarScan.v2.3.7.jar
   
scripts/create_baseline.sh:
>>> source /path/to/your/msings-virtual-environment/bin/activate
>>> VARSCAN=/path/to/your/msings-virtual-environment/bin/VarScan.v2.3.7.jar

scripts/run_msings.sh:
>>> source /path/to/your/msings-virtual-environment/bin/activate
>>> VARSCAN=/path/to/your/msings-virtual-environment/bin/VarScan.v2.3.7.jar

2. Run the create_intervals.sh bash script to create the msi_intervals file for your custom assay. This will create an msi_intervals file in the same directory as the bed file specified

>>> scripts/create_intervals.sh PATH/TO/BEDFILE

3. If necessary, bwa format and create a bwa index for your reference genome:

>>>  bwa index -a bwtsw ref_fasta

4. OPTIONAL: Now that we have CUSTOM_MSI_BED and CUSTOM_MSI_INTERVALS, you can update the create_baseline.sh script to point to these

>>> INTERVALS_FILE=/path/to/CUSTOM_MSI_INTERVALS;
>>> BEDFILE=/path/to/CUSTOM_MSI_BED;
>>> REF_GENOME=/path/to/REF_GENOME;

4. Create a file of the list of BAMS of MSI negative specimens, with each line being the absolute path to one sample

>>> /path/to/sampleA.bam
>>> /path/to/sampleB.bam
>>> /path/to/sampleC.bam

5. Create the baseline. Run the create_baseline.sh script for the batch of samples. Output will be in subfolders of the BAM data, subfolders named after the samples themselves

Default execution:
>>> scripts/create_baseline.sh PATH/TO/BAM_LIST PATH/TO/INTERVALS_FILE PATH/TO/BEDFILE PATH/TO/REF_GENOME

If you already edited the create_baseline.sh script to point to your reference files, you can instead just run:
>>> scripts/create_baseline.sh PATH/TO/BAM_LIST

RECOMMENDED: Now that the baseline file has been created, edit the baseline file to exclude loci which have standard deviations of zero. 
NOTE: Loci are excluded from the baseline file if the number of samples are insufficient to calculate statistics. 

The baseline contstruction process only need to be done once per assay/target data set. Files may be saved and re-used for subsequent analyses. 

9. OPTIONAL: Now we update the run_msings.sh to point to all the new custom files:

>>> INTERVALS_FILE=/path/to/CUSTOM_MSI_INTERVALS;
>>> BEDFILE=/path/to/CUSTOM_BEDFILE;
>>> MSI_BASELINE=/path/to/CUSTOM_MSI_BASELINE;
>>> REF_GENOME=/path/to/REF_GENOME;
 
10. Once the run_msings.sh script is updated for the new custom files, execution is the same as for Exome / TCGA data sets (above). 

Default execution:
>>>  scripts/run_msings.sh PATH/TO/BAM_LIST PATH/TO/INTERVALS_FILE PATH/TO/BEDFILE PATH/TO/REF_GENOME PATH/TO/MSI_BASELINE

If you already edited the create_baseline.sh script to point to your reference files:
>>>  scripts/run_msings.sh PATH/TO/BAM_LIST

 
Test to insure proper installation of scripts:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

>>>   cd msings
>>>   source msings-env/bin/active
>>>    ./testall
        Ran 11 tests in 0.068s
        OK

https://bitbucket.org/uwlabmed/msings
