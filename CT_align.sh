#!/bin/bash
#$ -pe smp 20
#$ -l h_vmem=2G
#$ -l h_rt=3:0:0
#$ -cwd 
#$ -j y
#$ -N Alignment_CT
#$ -m be

source /data/path/to/env/bin/activate

module load trimmomatic/
module load bowtie2/
module load samtools
module load bedtools

gunzip ./*.fq.gz
for file in ./*.fq
do
  base=`basename $file L1_1.fq`
  echo $file
  java -jar /data/Blizard-CGH-MadapuraLab/Debosree/CutnTag082020/fastqcinputs_082020/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  $file ${base}L1_2.fq ${base}trim_paired_L1_1.fq ${base}trim_unpaired_L1_1.fq ${base}trim_paired_L1_2.fq ${base}trim_unpaired_L1_2.fq \
  ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  bowtie2 -p 20 -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 \
  -x /data/Blizard-CGH-MadapuraLab/Debosree/CutnTag032020/hg38/hg38bowtie2 \
  -1 ${base}L1_1.fq -2 ./${base}L1_2.fq |
  samtools view -h -S -b -f 0x2 --threads 20 | 
  samtools sort -@ 20 > ${base}aln_sorted.bam 
  samtools index -b ${base}aln_sorted.bam
  samtools sort -@ 20 -n ${base}aln_sorted.bam > ${base}_byname.bam
  bedtools bamtobed -bedpe -i ${base}_byname.bam > ${base}.bed
  awk '$1==$4 && $6-$2 < 1000 {print $0}' ${base}.bed > ${base}.clean.bed
  cut -f 1,2,6 ${base}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${base}.clean.fragments.bed
  bedtools genomecov -bg -i ${base}.clean.fragments.bed -g ./hg38.chrom.sizes > ${base}.clean.fragments.bedgraph  
  bamCoverage -b ${base}aln_sorted.bam -o ${base}aln_sorted.bw --skipNAs --smoothLength 60 --centerReads -p 24 --normalizeUsing CPM
done