#!/bin/bash
#$ -pe smp 20
#$ -l h_vmem=2G
#$ -l h_rt=3:0:0
#$ -cwd 
#$ -j y
#$ -N Alignment_CT
#$ -m be

#need to install deeptools in your conda or python environment
#ml miniforge
#conda create -n Deeptools
#conda install -c bioconda deepTools
#conda acyivate Deeptools
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
  bedtools genomecov -bg -i ${base}.clean.fragments.bed -g ./hg38.chrom.sizes > ${base}.clean.fragments.bedgraph  #bedgraphs needed for peak calling using SEACR
  bamCoverage -b ${base}aln_sorted.bam -o ${base}aln_sorted.bw --skipNAs --smoothLength 60 --centerReads -p 24 --normalizeUsing CPM #bigwig files

  #Peak-calling using SEACR
  ml R
   bash path/to/SEACR/SEACR_1.3.sh ${base}.clean.fragments.bedgraph 0.000001 non relaxed ${base}_peaks #for relaxed peaks without IgG as background
      bash path/to/SEACR/SEACR_1.3.sh ${base}.clean.fragments.bedgraph 0.01 non stringent ${base}_peaks #for stringent peaks without IgG as background
done



------------------------------

#Plotting the data using deeptools

#!/bin/bash
#$ -pe smp 16
#$ -l h_vmem=2G
#$ -l h_rt=1:0:0
#$ -cwd
#$ -j y
#$ -m be

ml miniforge
conda activate Deeptools

computeMatrix reference-point --referencePoint TSS -S CUT_n_Tag1.bw CUT_n_Tag2.bw CUT_n_Tag_n.bw -bed Genes.bed -o Genes_CUT_n_Tag_sig.gz -a 5000 -b 5000 -bs 20 --skipZeros -p 20
#(-S bigwigs to use, -R bed files (regions) for plotting, -a nd -b distance from the TSS in bases, -p number of processors -bs binsize for counting signal)
plotHeatmap -m Genes_CUT_n_Tag_sig.gz -o Genes_CUT_n_Tag_sig.pdf --missingDataColor none -max value -min value --dpi 300 --colorList list_of_colors_to_plot or colorMap
#Generates heatmap from the matrix output of computeMatrix. 
plotProfile -m Genes_CUT_n_Tag_sig.gz -o Genes_CUT_n_Tag_sig_pro.pdf 
#Generates metasummary plot from the matrix output of computeMatrix. 

#QC_and_correlation_using_multiBam_or_multiBiwgwig_summary
multiBamSummary bins -b CUT_n_Tag1.bam CUT_n_Tag2.bam CUT_n_Tagn.bam -o CUT_n_Tag_mutliBam_summ.npz -p 16 --smartLabels
#generates correlation matrix for the bam ffiles used as input which can be plotted as PCA or correlation plots
 plotCorrelation --corData CUT_n_Tag_mutliBam_summ.npz --corMethod spearman or pearson --whatToPlot heatmap --skipZeros --plotNumbers -o CUT_n_Tag_mutliBam_summ_corr.pdf

 


 
