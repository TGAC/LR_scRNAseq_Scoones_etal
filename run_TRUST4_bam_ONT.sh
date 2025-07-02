#!/bin/bash
#SBATCH -p ei-medium # queue
#SBATCH -N 1 # number of nodes
#SBATCH -c 8 # number of cores
#SBATCH --mem 100G # memory pool for all cores
#SBATCH -J TRUST4_VDJ
#SBATCH -t 00-08:00 # time (D-HH:MM)
#SBATCH -o /ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/src/dump/TRUST4_ONT.%N.%j.out # STDOUT
#SBATCH -e /ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/src/dump/TURST4_ONT.%N.%j.err # STDERR

#source /jic/software/staging/RCSUPPORT-1775/stagingloader # TRUST4 testing install
source package 8ebd0f3f-543c-45f9-ac71-f8097550a9cb # TRUST4

# XM is UMI field in PacBio bam header
# UB is UMI field in 10X Illumina bam header

run-trust4 -t 8 --barcode CB --UMI XM -f /ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/references/hg38_bcrtcr.fa --ref /ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/references/human_IMGT+C.fa -b $1 --od $2 -o $3   
