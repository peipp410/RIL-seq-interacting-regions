#!/bin/bash

read1='./publications_data/ecoli/hfq_proq/hfq_lb/raw/hfq_lb_1.fastq.gz'
read2='./publications_data/ecoli/hfq_proq/hfq_lb/raw/hfq_lb_2.fastq.gz'
# Note! Please specify your own sequencing data.
gfasta='./ecoli_v38/ecoli_v38.fna'
gff='./annotation/GCF_000005845.2_ASM584v2_genomic.gff'

function rilseq {
   local read1=$1
   local read2=$2
   local gfasta=$3
   local gff=$4
   local prefix=$5

   if [ ! -d mapping ];then
      mkdir mapping
   fi

   map_single_fragments.py $gfasta -1 $read1 -2 $read2 -d mapping -r

   if [ ! -d readlist ];then
      mkdir readlist
   fi

   pre=`basename $read1 .gz`
   map_chimeric_fragments.py $gfasta mapping/${pre}_bwa.bam -a readlist/${prefix}.single_read.txt -A -r -d readlist >readlist/${prefix}.chimera_read.txt

   if [ ! -d result ];then
      mkdir result
   fi

   # if [ ! -d interaction ];then
   #    mkdir interaction
   # fi

   RILseq_significant_regions.py readlist/${prefix}.chimera_read.txt --servers 5 >result/${prefix}.inter.interaction.txt

   bwa mem -t 4 -R '@RG\tID:001\tPL:illumina\tSM:rep1' $gfasta $read1 | samtools sort -@ 4 -o mapping/${prefix}_R1.sorted.bam

   bwa mem -t 4 -R '@RG\tID:001\tPL:illumina\tSM:rep1' $gfasta $read2 | samtools sort -@ 4 -o mapping/${prefix}_R2.sorted.bam
   samtools index mapping/${prefix}_R1.sorted.bam
   samtools index mapping/${prefix}_R2.sorted.bam
 
   python3 annotation_motif.py -g $gfasta -t $gff -i ./result/${prefix}.inter.interaction.txt -c ./readlist/${prefix}.chimera_read.txt -b1 mapping/${prefix}_R1.sorted.bam -b2 mapping/${prefix}_R2.sorted.bam
   # see 'python3 annotation_motif.py -h' for detailed information.
}

rilseq $read1 $read2 $gfasta $gff hfq_lb

