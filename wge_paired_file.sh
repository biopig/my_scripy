#data=2018.5.18
#data1=2018.5.26
#user=zhanweimin
#630950832@qq.com
#在使用之前，现将代码阅读一遍，明确如何使用
#需要输入的内容，样品名称，每一对read只需要一个前缀就可以了，一个简单的例子：sh wgs.sh sample_name read1_profix ,假设文件的形式是  **.1_clean.fq.gz, **.2_clean.fq.gz
#sh wgs_paired_file.sh sample_name read1_profix read2_profix
#在使用时，可以根据实际情况修改线程数和内存大小

sample=$1
read1=$2
read2=$3
read3=$4
read4=$5
read5=$6
read6=$7

gatk=/home/zz/biosoft/gatk-4.0.4.0/gatk-package-4.0.4.0-local.jar
picard=/home/zz/miniconda3/share/picard-2.17.11-0/picard.jar
vep=/home/zz/biosoft/vep/ensembl-vep-release-92/vep

index=/media/zz/2545f6e0-b726-4a77-ac12-9942e94489a2/DNA-seq/index/Zea_mays.AGPv4_index
reference=/media/zz/2545f6e0-b726-4a77-ac12-9942e94489a2/DNA-seq/index/Zea_mays.AGPv4.dna.toplevel.fa

reference_dir=/media/zz/2545f6e0-b726-4a77-ac12-9942e94489a2/DNA-seq/index/

known_vcf=/media/zz/2545f6e0-b726-4a77-ac12-9942e94489a2/reference/zea_mays.vcf
known_vcf_dir=/media/zz/2545f6e0-b726-4a77-ac12-9942e94489a2/reference/

cache_dir=/media/zz/2545f6e0-b726-4a77-ac12-9942e94489a2/reference/vep/

if [ ! -d $sample/ ]
then mkdir ${sample}/ && cd ${sample}/
fi

#判断一些重要文件是否存在，不存在则创建文件
if [ ! -e ${reference_dir}/Zea_mays.AGPv4.dna.toplevel.dict ]
then java -jar -Xmx5G $picard CreateSequenceDictionary R=$reference O=${reference_dir}/Zea_mays.AGPv4.dna.toplevel.dict  #为参考序列生成一个dict文件
fi

if [ ! -e ${known_vcf_dir}/zea_mays.vcf.idx ]
then java -jar -Xmx5G $gatk IndexFeatureFile -F $known_vcf  #为参考的变异文件构建索引
fi

if [ ! -e ${reference_dir}/Zea_mays.AGPv4_index.1.bt2 ]
then bowtie2-build ${reference} Zea_mays.AGPv4_index 1>Zea_mays.AGPv4.bowtie2_index.log 2>&1 #构建索引
fi

#判断你的文件有几组，然后进行比对
if [ $#==2 ]
then bowtie2 -p 4 -x $index -1 ${read1}_1.fastq.gz -2 ${read1}_2.fastq.gz -S ${sample}.sam 1>bowtie2_alignment.log 2>bowtie_alignment.err
elif [ $#==3 ]
then bowtie2 -p 4 -x $index -1 ${read1}_1.fastq.gz -2 ${read1}_2.fastq.gz -1 ${read2}_1.fastq.gz -2 ${read2}_2.fastq.gz -S ${sample}.sam 1>bowtie2_alignment.log 2>bowtie_alignment.err
elif [ $#==4 ]
then bowtie2 -p 4 -x $index -1 ${read1}_1.fastq.gz -2 ${read1}_2.fastq.gz -1 ${read2}_1.fastq.gz -2 ${read2}_2.fastq.gz -1 ${read3}_1.fastq.gz -2 ${read3}_2.fastq.gz -S ${sample}.sam 1>bowtie2_alignment.log 2>bowtie_alignment.err
elif [ $#==5 ]
then bowtie2 -p 4 -x $index -1 ${read1}_1.fastq.gz -2 ${read1}_2.fastq.gz -1 ${read2}_1.fastq.gz -2 ${read2}_2.fastq.gz -1 ${read3}_1.fastq.gz -2 ${read3}_2.fastq.gz -1 ${read4}_1.fastq.gz -2 ${read4}_2.fastq.gz -S ${sample}.sam 1>bowtie2_alignment.log 2>bowtie_alignment.err
elif [ $#==6 ]
then bowtie2 -p 4 -x $index -1 ${read1}_1.fastq.gz -2 ${read1}_2.fastq.gz -1 ${read2}_1.fastq.gz -2 ${read2}_2.fastq.gz -1 ${read3}_1.fastq.gz -2 ${read3}_2.fastq.gz -1 ${read4}_1.fastq.gz -2 ${read4}_2.fastq.gz -1 ${read5}_1.fastq.gz -2 ${read5}_2.fastq.gz -S ${sample}.sam 1>bowtie2_alignment.log 2>bowtie_alignment.err
elif [ $#==7 ]
then bowtie2 -p 4 -x $index -1 ${read1}_1.fastq.gz -2 ${read1}_2.fastq.gz -1 ${read2}_1.fastq.gz -2 ${read2}_2.fastq.gz -1 ${read3}_1.fastq.gz -2 ${read3}_2.fastq.gz -1 ${read4}_1.fastq.gz -2 ${read4}_2.fastq.gz -1 ${read5}_1.fastq.gz -2 ${read5}_2.fastq.gz -1 ${read6}_1.fastq.gz -2 ${read6}_2.fastq.gz -S ${sample}.sam 1>bowtie2_alignment.log 2>bowtie_alignment.err
fi

#对sam文件进行排序
samtools view -bS -@ 4 ${sample}.sam |samtools sort -@ 4 - -T ${sample}.sorted -o ${sample}.sorted.bam &&\

#对bam文件构建索引
samtools index ${sample}.sorted.bam &&\

#标记重复
java -jar -Xmx5G $gatk MarkDuplicates -I ${sample}.sorted.bam -M ${sample}.markdup_matrics.txt -O ${sample}.sorted.markdup.bam &&\

#对bam文件加头文件
java -jar -Xmx5G $picard AddOrReplaceReadGroups I=${sample}.sorted.markdup.bam O=${sample}.addhead.bam RGLB=${sample}ID RGPL=illumina RGPU=${sample}PU RGSM=${sample} &&\

#碱基质量的重矫正
java -jar -Xmx5G $gatk BaseRecalibrator -I ${sample}.addhead.bam --known-sites $known_vcf -O ${sample}_baserecalibrator.list -R $reference && \
java -jar -Xmx5G $gatk ApplyBQSR -bqsr ${sample}_baserecalibrator.list -I ${sample}.addhead.bam -O ${sample}_applybqsr.bam && \

#变异的检测，生成vcf文件
java -jar -Xmx5G $gatk HaplotypeCaller -I ${sample}_applybqsr.bam -O ${sample}_haplotypecaller.vcf -R $reference &&\

#这一步就是对变异位点进行质量的重矫正，可是进行下面这两步后的变异数目没有变化，需要继续探究
#java -jar -Xmx4G $gatk VariantRecalibrator -R $reference -V ${sample}.haplotypecaller.vcf -resource hapmap,known=false,training=true,truth=true,prior=15.0:$known_vcf -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 --tranches-file ${sample}.variantrecalibrator.snps.tranches -O ${sample}.snps.recal --max-gaussians 4
#java -jar -Xmx4G $gatk ApplyVQSR -R $reference -V ${sample}.haplotypecaller.vcf --ts-filter-level 99.0 --tranches-file ${sample}.variantrecalibrator.snps.tranches --recal-file ${sample}.snps.recal -mode SNP -O ${sample}.snps.VQSR.vcf

#对变异文件vcf进行注释
$vep --fasta $reference --offline --species zea_mays --cache_version 39 --dir_cache $cache_dir -i ${sample}_haplotypecaller.vcf -o ${sample}.vep.annotate.vcf


