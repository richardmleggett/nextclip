#!/bin/bash

export ref_min_size=0
export libdir=examples/LIB1468
export lib=LIB1468
export read_one=LIB1468part_ATCACG_L001_R1_001.fastq
export read_two=LIB1468part_ATCACG_L001_R2_001.fastq
export reference=examples/Reference/AL645882.fasta
export organism="Streptomyces coelicolor"
scripts/nextclip_lmp_analysis.sh
