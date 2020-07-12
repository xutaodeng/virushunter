#!/bin/bash

sample=$1
thread=12

STAR=/mnt/cluster/xdeng/tools/STAR-STAR_2.4.0k/bin/Linux_x86_64/STAR
genomeDir=${PWD}/hg19_1/
ref=/mnt/cluster/xdeng/gatk2.8Resource/human_g1k_v37.fasta
runDir1=${PWD}/raw_data/${sample}/
runDir2=${PWD}/raw_data/${sample}/2pass/
genomeDir2=${PWD}/raw_data/${sample}/hg19_2/
read1=${PWD}/raw_data/fastq/${sample}_R1.fastq.gz
read2=${PWD}/raw_data/fastq/${sample}_R2.fastq.gz

########################################
# STAR
# #########################################
#only run the first time to generate genome index
#echo $STAR --runMode genomeGenerate --genomeDir $genomeDir --genomeFastaFiles $ref  --runThreadN ${thread}
echo mkdir $runDir1
echo cd $runDir1
echo $STAR --genomeDir $genomeDir --readFilesIn $read1 $read2 --readFilesCommand zcat --runThreadN ${thread}
echo mkdir $genomeDir2
echo $STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $ref --sjdbFileChrStartEnd ${runDir1}SJ.out.tab --sjdbOverhang 75 --runThreadN ${thread}
echo mkdir $runDir2
echo cd $runDir2
echo $STAR --genomeDir $genomeDir2 --readFilesIn $read1 $read2 --readFilesCommand zcat --runThreadN ${thread}

########################################
# GATK
#########################################
s1=${runDir2}Aligned.out.sam
s2=${runDir2}s2.bam
s3=${runDir2}s3.bam
s4=${runDir2}s4.bam
s5=${runDir2}s5.bam
s6=${runDir2}s6.bam
s7=${runDir2}s7.bam
out1=${runDir2}output1.vcf
out2=${runDir2}output2.vcf
out3=${runDir2}output3.vcf
outhtml=${runDir2}output3.html
recal=${runDir2}recal.table
dupfile=${runDir2}output.metrics

picard=/mnt/cluster/xdeng/tools/picard-tools-1.119/
gatk=/mnt/cluster/xdeng/tools/GenomeAnalysisTK-3.3-0/
dbSNP=/mnt/cluster/xdeng/gatk2.8Resource/dbsnp_138.b37.vcf
knownindel=/mnt/cluster/xdeng/gatk2.8Resource/Mills_and_1000G_gold_standard.indels.b37.vcf

option='-Djava.io.tmpdir=/mnt/cluster2/tmp  -Xms1024m -Xmx4096m '
echo java $option -jar ${picard}SortSam.jar I=${s1} O=${s2} SO=coordinate VALIDATION_STRINGENCY=SILENT
echo java  $option -jar ${picard}AddOrReplaceReadGroups.jar I=${s2} O=${s3} SO=coordinate RGID=id RGLB=solexa-123 RGPL=ILLUMINA RGPU=AXL2342 RGSM=${sample}
echo java $option -jar ${picard}MarkDuplicates.jar I=${s3} O=${s4}  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${dupfile}
echo samtools index $s4
echo java $option -jar ${gatk}GenomeAnalysisTK.jar -T SplitNCigarReads -R $ref -I ${s4} -o ${s5} -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
echo java $option -jar ${gatk}GenomeAnalysisTK.jar -T BaseRecalibrator -nct ${thread}  -R $ref  -I $s5 -knownSites $dbSNP -knownSites $knownindel  -o $recal
echo java $option -jar ${gatk}GenomeAnalysisTK.jar -T PrintReads -nct ${thread} -R $ref -I $s5 -BQSR $recal -o $s6
echo java $option -jar ${gatk}GenomeAnalysisTK.jar -T HaplotypeCaller -R $ref -I $s6 --dbsnp $dbSNP -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o $out1
echo java $option -jar ${gatk}GenomeAnalysisTK.jar -T VariantFiltration -R $ref -V $out1 -window 35 -cluster 3 -filterName FS -filter \"FS ">" 30.0\" -filterName QD -filter \"QD "<" 2.0\" -o $out2
echo java $option -jar  /mnt/cluster/tools/snpEff/snpEff.jar eff -c /mnt/cluster/tools/snpEff/snpEff.config -s outhtml -v GRCh37.74 $out2 ">" $out3