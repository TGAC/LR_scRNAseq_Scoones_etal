#!/bin/bash
#SBATCH -p ei-medium # queue
#SBATCH -N 1 # number of nodes
#SBATCH -c 20 # number of cores
#SBATCH -J samstats # job name
#SBATCH -t 2-00:00 # time (D-HH:MM)
#SBATCH --mem 100G # memory pool for all cores
#SBATCH -o /ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/logs/samstats.%N.%j.out # STDOUT
#SBATCH -e /ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/logs/samstats.%N.%j.err # STDERR

source package c92263ec-95e5-43eb-a527-8f1496d56f1a #samtools-1.18

# inbam=$1
# class_file=$2
# workdir=$(dirname $class_file)

inbam="/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/data/MASseq/nanopore/wf-single-cell-YL/Cribbs.run2/Earlham1BPBMC_Ian/bams/Earlham1PBMC_Iain.merged.sorted.bam"
class_file="/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/data/MASseq/nanopore/wf-single-cell-YL/Cribbs.run2/Earlham1BPBMC_Ian/sqanti3_output/default_filter/default_filter_RulesFilter_result_classification.txt"
read_id_file="/ei/projects/f/f4e517e2-faa5-488e-82ef-66c09af5aff9/data/data/MASseq/nanopore/wf-single-cell-YL/Cribbs.run2/Earlham1BPBMC_Ian/isoquant_output/OUT/OUT.transcript_model_reads.tsv.gz"
workdir=$(dirname $class_file)
cd $workdir

samtools flagstat -@ 20 $inbam > ${inbam%.bam}.flagstat.txt

echo -e "Class\tTotal\tArtifact" > sqanti3_classes.txt
awk 'NR > 1 {print $6}' $class_file | sort | uniq | while read value; do 
    total=$(grep -w "$value" $class_file | wc -l);
    artifact=$(awk -v val="$value" '$6 == val && $NF == "Artifact"' $class_file | wc -l);
    echo -e "$value\t$total\t$artifact";
done >> sqanti3_classes.txt

source /hpc-home/kudashev/.bashrc
conda activate lr_analysis #env with pandas 

python /ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/scripts/merge_sqanti_readid.py \
--class_file $class_file --read_id_file $read_id_file \
--output $workdir/read_id2sqanti_cat.tsv --filtered_readid $workdir/filtered_readid.txt

samtools view -@ 20 -b -N $workdir/filtered_readid.txt $inbam > ${inbam%.bam}.filtered.bam

samtools flagstat -@ 20 ${inbam%.bam}.filtered.bam > ${inbam%.bam}.filtered.flagstat.txt
