#!/bin/bash

mkdir demult2

sbatch --time 12:00:00 --mem 32G --wrap "/Users/simeon01/applications/demultiplexer.rhel/demuxFQ -c -d -i -e -t 1 -r 0.01 -R -l 9 -o ./demult2 -b ./demult2/SLX-21389.lostreads.r_1.fq.gz -s SLX-21389_summary_demuxFQ_r1.txt SLX-21389_demuxFQ_r1.txt SLX-21389.UnspecifiedIndex.H7J5VBGXL.s_1.r_1.fq.gz"

sbatch --mem 32G --time 12:00:00 --wrap "/Users/simeon01/applications/demultiplexer.rhel/demuxFQ -c -d -i -e -t 1 -r 0.01 -R -l 9 -o ./demult2 -b ./demult2/SLX-21389.lostreads.r_2.fq.gz -s SLX-21389_summary_demuxFQ_r2.txt SLX-21389_demuxFQ_r2.txt SLX-21389.UnspecifiedIndex.H7J5VBGXL.s_1.r_2.fq.gz"





#!/bin/bash
for f in *r1.fq.gz
do
  sbatch --time 01:00:00 --mem 5G --wrap "zcat $f | split -dl 4000000 - ${f%%fq.gz} --a 3 --additional-suffix=.fastq && zcat ${f/.r1./.r2.} | split -dl 4000000 - ${f/r1.fq.gz/r2.} --a 3 --additional-suffix=.fastq"
done



#!/bin/bash
mkdir trimmed2
# paired-end cut and tag
for file in *.r1.*.fastq
do
fq1=$file
fq2=${fq1/r1/r2}
ls $fq1
ls $fq2
echo "----"
sbatch --time 01:00:00 -o %j.out -e %j.err --mem 16000 --wrap "cutadapt -q 20  -o trimmed2/${fq1%%.fastq}.trimmed.fq.gz -p trimmed2/${fq2%%.fastq}.trimmed.fq.gz $fq1 $fq2 "
done




#!/bin/bash
mkdir aligned
g='~/hg38_ecoliK12/hg38_selected_ecoli_K12.fa'
w='~/hg38.whitelist.sorted.bed'

#aln with bwa
for file in *.r1*trimmed.fq.gz
  do
  f1=$file
  f2=${file/r1/r2}
  sbatch --time 12:00:00 --mem 16G --wrap "bwa mem -M $g $f1 $f2  | samtools view -Sb -F780 -q 10 -L $w - | samtools sort -@ 8 -   > aligned/${f1%%.trimmed.fq.gz}.hg38.sort.bam"
done


#!/bin/bash

for file in *.r1.00.hg38.sort.bam
do
sample=${file%%.r1.00.hg38.sort.bam}
sbatch --time 01:00:00  --mem 8G --wrap "samtools merge $sample.merged.bam ${file%%.r1.00.hg38.sort.bam}.r1.*.hg38.sort.bam"
echo "-------"
done
#!/bin/bash

for file in *.r1.00.hg38.sort.bam
do
sample=${file%%.r1.00.hg38.sort.bam}
sbatch --time 01:00:00  --mem 8G --wrap "samtools merge $sample.dupl.bam ${file%%.r1.00.hg38.sort.bam}.r1.*.hg38.sort.bam"
echo "-------"
done


#!/bin/bash

for bam in *merged.bam
do
sbatch --time 12:00:00 --mem 8G --wrap "java -Xmx7g -jar /home/simeon01/applications/picard-2.20.3.jar MarkDuplicates INPUT=$bam OUTPUT=${bam%%.bam}.markduplicates.bam REMOVE_DUPLICATES=true  AS=true METRICS_FILE=${bam%%.bam}.markduplicates_metrics.txt | samtools index ${bam%%.bam}.markduplicates.bam"
done

#!/bin/bash


for file in *merged.bam
do
bam_hg38=$file
bam_hg38_nodup=${file%%.bam}.markduplicates.bam

sbatch --time 00:05:00 --mem 8G --wrap "samtools view -c -F 260 $bam_hg38 >  ${file%%merged.bam}dupl.stat2"
sbatch --time 00:05:00 --mem 8G --wrap "samtools view -c -F 260 $bam_hg38_nodup >  ${file%%.bam}.stat5"
done

# collect stats by chr
for file in *.merged*bam
do
sbatch --mem 4G --time 00:05:00 --wrap  "samtools flagstat $file -O tsv > ${file%%.bam}.flagstat"
sbatch --mem 4G --time 00:05:00 --wrap  "samtools idxstats $file > ${file%%.bam}.idxstats"
done

#!/bin/bash

for bam in *.bam; do sbatch --time 01:00:00 --mem 4G --wrap "samtools index $bam"; done


# generate bw
 for bam in *.markduplicates.bam;
 do
 tot_r_hg38=`cat ${bam%%.markduplicates.bam}.stat5`;
 scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }');
 echo $scal_factor_hg38;
 echo sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw";
 sbatch --time 06:00:00  -e bw_gen%j.${bam%%.markduplicates.bam}.out --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.markduplicates.bam}.w10.rpm.bw"
 done

for bam in *.dupl.bam;
 do
 tot_r_hg38=`cat ${bam%%.bam}.stat2`;
 scal_factor_hg38=$(awk -v m=$tot_r_hg38 'BEGIN { print 1000000/m }');
 echo $scal_factor_hg38;
 echo sbatch --time 06:00:00 --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.bam}.w10.rpm.bw";
 sbatch --time 06:00:00  -e bw_gen%j.${bam%%.markduplicates.bam}.out --mem 16G --wrap "bamCoverage --scaleFactor $scal_factor_hg38 -bs 10 -b $bam -of \"bigwig\" -o ${bam%%.bam}.w10.rpm.bw"
 done



for file in *dupl.bam
do
sbatch --time 01:00:00 --wrap "bamCoverage --bam $file -o ${file%%bam}.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads"

done

for file in *dupl.bam
do
sbatch --time 01:00:00 --wrap "bamCoverage --bam $file -o ${file%%bam}.rpkm.bs10.bw --binSize 10 --normalizeUsing RPKM --extendReads"

done


for file in *.markduplicates.bam
do
sbatch --time 01:00:00 --wrap "bamCoverage --bam $file -o ${file%%bam}.cpm.bs10.bw --binSize 10 --normalizeUsing CPM  --extendReads"

done

for file in *.markduplicates.bam
do
sbatch --time 01:00:00 --wrap "bamCoverage --bam $file -o ${file%%bam}.rpkm.bs10.bw --binSize 10 --normalizeUsing RPKM --extendReads"

done

#!/bin/bash


# sort by name
# convert to bedpe
for file in *.merged.markduplicates.bam
do
label=${file%%.merged.markduplicates.bam}
cmd1="samtools sort -n $file > $label.sortName.bam &&\
  bedtools bamtobed -bedpe -i $label.sortName.bam > $label.sortName.bed"
echo $cmd1
sbatch --time 01:00:00  --mem 4G --wrap "$cmd1"
done

#!/bin/bash

# fragment size selection
genome=~/hg38_selected.sorted.genome
for file in *sortName.bed
do
cmd_1000="awk '{if (\$1==\$4 && \$6-\$2 < 1000) print \$0}' $file > ${file%%.bed}.1000.clean.bed && cut -f 1,2,6 ${file%%.bed}.1000.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${file%%.bed}.1000.clean.fragments.bed && bedtools genomecov -bg -i ${file%%.bed}.1000.clean.fragments.bed -g $genome > ${file%%.bed}.1000.clean.fragments.bedgraph"
echo $cmd_1000
echo " >> =============== ============= << "
sbatch --time 01:00:00 --mem 1G --wrap "$cmd_1000"
done

mkdir seacr_no_ctrl

for file in *sortName.bed
do
bdg_1000=${file%%.bed}.1000.clean.fragments.bedgraph
cmd_s_1000="~/applications/SEACR_1.3.sh $bdg_1000 0.05 non stringent seacr_no_ctrl/${bdg_1000%%clean}.0.05fdr"
echo $cmd_s_1000
sbatch --mem 1G --time 01:00:00 --wrap "$cmd_s_1000"
done


for file in *sortName.bed
do
bdg_1000=${file%%.bed}.1000.clean.fragments.bedgraph
cmd_s_1000="~/applications/SEACR_1.3.sh $bdg_1000 0.05 non relaxed seacr_no_ctrl/${bdg_1000%%clean}.rel.0.05fdr"
echo $cmd_s_1000
sbatch --mem 1G --time 01:00:00 --wrap "$cmd_s_1000"
done

#!/bin/bash

mkdir seacr_no_ctrl_01

for file in *sortName.bed
do
bdg_1000=${file%%.bed}.1000.clean.fragments.bedgraph
cmd_s_1000="~/applications/SEACR_1.3.sh $bdg_1000 0.01 non stringent seacr_no_ctrl_01/${bdg_1000%%clean}.0.01fdr"
echo $cmd_s_1000
sbatch --mem 1G --time 01:00:00 --wrap "$cmd_s_1000"
done


for file in *sortName.bed
do
bdg_1000=${file%%.bed}.1000.clean.fragments.bedgraph
cmd_s_1000="~/applications/SEACR_1.3.sh $bdg_1000 0.01 non relaxed seacr_no_ctrl_01/${bdg_1000%%clean}.rel.0.01fdr"
echo $cmd_s_1000
sbatch --mem 1G --time 01:00:00 --wrap "$cmd_s_1000"
done

#overalps of technical replicates
#!/bin/bash

frag_cases=(DOX3_200nM DOX3_20nM DOX3_2nM )
B=( 1 2 )
T=(1 2 3 4 5 )
min_reads=(0 5 8 10)

for file in ${frag_cases[@]}
do
for b in ${B[@]}
do
for read in ${min_reads[@]}
do

intervene venn --save-overlaps -i *${file}*B${b}*min${read}.bed -o Intervenn_${file}_B${b}_min${read}
intervene upset  --save-overlaps -i *${file}*B${b}*min${read}.bed -o Intervenn_${file}_B${b}_min${read}
intervene pairwise  --save-overlaps -i *${file}*B${b}*min${read}.bed -o Intervenn_${file}_B${b}_min${read}

multiIntersectBed -wa -wb  -i *${file}*B${b}*min${read}.bed | awk '{if ($4>=1) print $0}' > ${file}.B${b}.min${read}.multi1.bed
multiIntersectBed -wa -wb  -i *${file}*B${b}*min${read}.bed | awk '{if ($4>=2) print $0}' > ${file}.B${b}.min${read}.multi2.bed
multiIntersectBed -wa -wb  -i *${file}*B${b}*min${read}.bed | awk '{if ($4>=3) print $0}' > ${file}.B${b}.min${read}.multi3.bed
multiIntersectBed -wa -wb  -i *${file}*B${b}*min${read}.bed | awk '{if ($4>=4) print $0}' > ${file}.B${b}.min${read}.multi4.bed
multiIntersectBed -wa -wb  -i *${file}*B${b}*min${read}.bed | awk '{if ($4>=5) print $0}' > ${file}.B${b}.min${read}.multi5.bed

done
done
done

#merge bed files
 for f in *bed ; do mergeBed -i $f > ${f%%.bed}.m.bed ; done

 for file in ${frag_cases[@]}
 do

 min_reads=(0 5 8 10)
 multi=(1 2 3 4 5 )
 for read in ${min_reads[@]}
 do
 for m in ${multi[@]}
 do
 intervene venn  -i  *.min${read}.multi${m}.bed -o Intervenn_${file}_Ba_min${read}_multi${m}
 intervene pairwise  -i  *.min${read}.multi${m}.bed -o Intervenn_${file}_Ba_min${read}_multi${m}
 intervene upset -i *.min${read}.multi${m}.bed -o Intervenn_${file}_Ba_min${read}_multi${m}

 done
 done
 done


 sbatch --time 2-12:00:00 --mem 8G --ntasks 20 --wrap "multiBigwigSummary BED-FILE --BED 11_11_concensusbed.bed  -b mol1*mol2*bw -o results_mol1_mol2_BED.npz  -p 20 "

 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b1_BED.npz --corMethod pearson --whatToPlot heatmap --skipZeros --colorMap RdYlBu --plotNumbers -o Pearson_heat_nz.BRD4_B1_BED.pdf"
 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b2_BED.npz --corMethod spearman --whatToPlot heatmap --skipZeros --colorMap RdYlBu --plotNumbers -o Spearman_heat_nz.BRD4.B1_BED.pdf"
 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b1_BED.npz --corMethod pearson --whatToPlot heatmap --skipZeros --colorMap RdYlBu --plotNumbers -o Pearson_heat_nz.BRD4_B2_BED.pdf"
 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b2_BED.npz --corMethod spearman --whatToPlot heatmap --skipZeros --colorMap RdYlBu --plotNumbers -o Spearman_heat_nz.BRD4.B2_BED.pdf"


 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b1_BED.npz --corMethod pearson  --skipZeros --whatToPlot scatterplot  -o Pearson_scatterBRD4_nz_B1_BED.pdf  "
 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b2_BED.npz --corMethod spearman --skipZeros --whatToPlot scatterplot  -o Spearman_scatterBRD4_nz_B1_BED.pdf  "
 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b1_BED.npz --corMethod pearson  --skipZeros --whatToPlot scatterplot  -o Pearson_scatterBRD4_nz_B2_BED.pdf  "
 sbatch --time 12:00:00 --mem 12G --ntasks 10  --wrap "plotCorrelation --corData results_BRD4_b2_BED.npz --corMethod spearman --skipZeros --whatToPlot scatterplot  -o Spearman_scatterBRD4_nz_B2_BED.pdf  "
