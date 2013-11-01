Examples

This directory contains some example data for testing NextClip. This consists of:
1. Inside the reads directory: 125,000 pairs of reads from a (not especially high quality) Streptomyces coelicolor Nextera LMP run.
2. Inside the reference directory: the S. coelicolor reference, from https://www.ebi.ac.uk/ena/data/view/AL645882.
3. Inside the configure directory, LIB468.txt is a parameters file for passing to the NextClip pipeline.

Running NextClip standalone

First compile NextClip as detailed in the manual. To run the tool, change into the directory containing the reads, configure, reference etc. directories and type:

nextclip -i reads/LIB1468part_ATCACG_L001_R1_001.fastq -j reads/LIB1468part_ATCACG_L001_R2_001.fastq -n 125000 -o output

Running the NextClip analysis pipeline

To run the NextClip pipeline, ensure you have BWA, R and LaTeX installed and available and that you have carried out all the other installation steps detailed in the manual.

scripts/nextclip_lmp_analysis.pl -config configure/LIB1468.txt -scheduler none -bwathreads 1

If all works successfully, a new directory called LIB1468 should appear and inside that will be a latex directory containing a PDF file.