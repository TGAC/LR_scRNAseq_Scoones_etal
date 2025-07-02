#!/bin/bash
#SBATCH -p ei-long # queue
#SBATCH -N 1 # number of nodes
#SBATCH -c 20 # number of cores
#SBATCH -J MASseq_SQANTI3 # job name
#SBATCH --mem 250G # memory pool for all cores
#SBATCH -t 8-00:00 # time (D-HH:MM)
#SBATCH -o /ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/logs/MASseq_SQANTI3.%N.%j.out # STDOUT
#SBATCH -e /ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/logs/MASseq_SQANTI3.%N.%j.err # STDERR

## DEFINE GLOBAL VARIABLES
workdir="/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/data/MASseq/nanopore/wf-single-cell-YL/Cribbs.run2/Earlham1BPBMC_Ian"
ref_gtf="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/refdata-gex-GRCh38-2024-A/genes/genes.gtf"
ref_fa="/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/references/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
bam="/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/data/MASseq/nanopore/wf-single-cell-YL/Cribbs.run2/Earlham1BPBMC_Ian/bams/Earlham1PBMC_Iain.merged.sorted.bam"

isoquant="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/tools/isoquant"

singularity exec $isoquant/isoquant.img $isoquant/IsoQuant/isoquant.py  -t 20 --data_type nanopore \
    --reference $ref_fa --genedb $ref_gtf --bam $bam --complete_genedb \
    --splice_correction_strategy conservative_ont \
    --matching_strategy  default \
    --model_construction_strategy sensitive_ont \
    --polya_requirement auto --report_canonical all \
    --output $workdir/isoquant_output --sqanti_output

isoquant_gtf="$workdir/isoquant_output/OUT/OUT.transcript_models.strand_fixed.gtf"

source /hpc-home/kudashev/.bashrc
conda activate SQANTI3.env # SQANTI3-5.2.1
SQANTI3="/ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/tools/SQANTI3-5.2.1"

singularity exec $SQANTI3/SQANTI3-5.2.1.img sqanti3_qc.py --cpus 20 --report pdf --force_id_ignore --aligner_choice minimap2 \
--polyA_motif_list $SQANTI3/data/polyA_motifs/mouse_and_human.polyA_motif.txt \
--CAGE_peak $SQANTI3/data/ref_TSS_annotation/human.refTSS_v3.1.hg38.bed \
--sites GTAG,GCAG,ATAC --dir $workdir/sqanti3_output/singularity_test --ratio_TSS_metric max \
--SR_bam /ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/data/MASseq/illumina/Data_Package_Batch_2023_01_18_illumina/CB-GENANNO-551_Iain_Macaulay_EI_AL_ENQ-5477_A_01/bamfile_Iain_Macaulay_EI_AL_ENQ-5477_A_01_Analysis/bam_file_paths.txt \
$isoquant_gtf $ref_gtf $ref_fa

singularity exec $SQANTI3/SQANTI3-5.2.1.img sqanti3_filter.py rules --skip_report -j $workdir/sqanti3_output/filter_rules.json \
--gtf $workdir/sqanti3_output/OUT.transcript_models.strand_fixed_corrected.gtf $workdir/sqanti3_output/OUT.transcript_models.strand_fixed_classification.txt 

singularity exec $SQANTI3/SQANTI3-5.2.1.img sqanti3_filter.py rules --skip_report -j $workdir/sqanti3_output/default_rules.json --output default_filter \
--gtf $workdir/sqanti3_output/OUT.transcript_models.strand_fixed_corrected.gtf $workdir/sqanti3_output/OUT.transcript_models.strand_fixed_classification.txt 

python /ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/scripts/merge_sqanti_readid.py \
--class_file $workdir/sqanti3_output/OUT.transcript_models.strand_fixed_classification.txt --read_id_file $workdir/isoquant_output/OUT/OUT.transcript_model_reads.tsv.gz \
--output $workdir/read_id2sqanti_cat.tsv
