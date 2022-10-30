
#bamcompare
file=KATO3
samtools index ./$file/$file'_H3K27ac_picard_.sort.bam'
samtools index ./$file/$file'_Input_picard_.sort.bam'
bamCompare -p 10 -b1 ./$file/$file'_H3K27ac_picard_.sort.bam' -b2 ./$file/$file'_Input_picard_.sort.bam' --skipNAs --scaleFactorsMethod readCount --operation subtract --outFileFormat bedgraph -o ./$file/$file'.bedgraph' --extendReads 200

awk '{if($4<0){$4=0};print}' ./$file/$file'.bedgraph' > ./$file/$file'_move0.bedgraph'
totalreads= awk '{sum=sum+$4}END{print sum}' ./$file/$file'_move0.bedgraph'
awk -v totalreads=6.97793e+07 '{$4=$4/totalreads*1000000;print}' ./$file/$file'_move0.bedgraph' > ./$file/$file'_rpm.bedgraph'
sort -k1,1 -k2,2n ./$file/$file'_rpm.bedgraph' > ./$file/$file'_sort.bedgraph'

bedGraphToBigWig ./$file/$file'_sort.bedgraph' hg19.fa.fai ./$file/$file'.bw'


####bamcoverage
bamCoverage --bam N_out_all.bam -o N_out_all.bw  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize 2864785220 --ignoreForNormalization chrX --extendReads	

###

computeMatrix reference-point \
 -S TP63.sort.bam.ym.bw \
 -R promoter.bed enhancer.bed \
 --referencePoint center \
 -a 2000 -b 2000 -out TP63.gz
 
 

 plotHeatmap -m TP63.gz \
     -out TP63.pdf \
     --colorMap Oranges Blues BrBG copper_r \
     --whatToShow 'plot,heatmap and colorbar' 