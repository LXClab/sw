#conda activate chip
#do
#双端
file=Fadu
cd /home/jyxyyxy/publicNAS1/YM-stomach/
#input组
bowtie2 -p 100 -x ./bowtie2index/hg19/hg19 -1 ./Fadu/Fadu_In_clean_1.fq.gz -2 ./Fadu/Fadu_In_clean_2.fq.gz -S ./Fadu/$file'_Input.sam'
samtools view -@ 10 -b -S ./$file/$file'_Input.sam' > ./$file/$file'_Input.bam'
samtools sort -@ 10 ./$file/$file'_Input.bam' -o ./$file/$file'_Input.sort.bam'
samtools index ./$file/$file'_Input.sort.bam'

#H3K27ac组或者IP组
bowtie2 -p 100 -x /home/jyxyyxy/bowtie2index/hg19/hg19 -1 ./Fadu/Fadu_IP_clean_1.fq.gz -2 ./Fadu/Fadu_IP_clean_2.fq.gz -S ./Fadu/$file'_IP.sam'
samtools view -@ 10 -b -S ./$file/$file'_IP.sam' > ./$file/$file'_IP.bam'
samtools sort -@ 10 ./$file/$file'_IP.bam' -o ./$file/$file'_IP.sort.bam'
samtools index ./$file/$file'_IP.sort.bam'

##去除PCR重复
picard MarkDuplicates \INPUT=./$file/$file'_Input.sort.bam' \OUTPUT=./$file/$file'_Input_picard_.sort.bam' \METRICS_FILE= ./$file/$file'.Input_markDup1.metric' REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 
picard MarkDuplicates \INPUT=./$file/$file'_IP.sort.bam' \OUTPUT=./$file/$file'_Ip_picard_.sort.bam' \METRICS_FILE= ./$file/$file'.Ip_markDup1.metric' REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=true 


#call peak
macs2 callpeak  -t ./$file/$file'_Ip_picard_.sort.bam' -c ./$file/$file'_Input_picard_.sort.bam' -f BAMPE -g hs --keep-dup all --name ./$file/$file --bdg --broad --SPMR --nomodel --extsize 200 -q 0.01
macs2 bdgcmp -t ./$file/$file'_treat_pileup.bdg' -c ./$file/$file'_control_lambda.bdg' -o ./$file/$file'.bdg' -m FE


awk '{print $1,$2,$3,$4,$5}' ./$file/$file'_peaks.broadPeak' > ./$file/$file'_broadPeak.bed'
sed -i 's/ /\t/g' ./$file/$file'_broadPeak.bed'
awk 'BEGIN{OFS=" "}{$5= "@"}{$6= $4}{$7 = "."}{$8 = "@"}{$9 = "@"}{print $1,$4,$5,$2,$3,$8,$7,$9,$6}' ./$file/$file'_broadPeak.bed' > ./$file/$file'_broadPeak.gff'
sed -i 's/@//g' ./$file/$file'_broadPeak.gff'
sed -i 's/ /\t/g' ./$file/$file'_broadPeak.gff'

mkdir ./$file/rose-results


