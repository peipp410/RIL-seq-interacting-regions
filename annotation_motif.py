import argparse
import pysam
import pandas as pd
from collections import defaultdict
import os


def read_gtf(filename, keep_feature=None, name='gene'):
    """
    Read in an annotation file (gff or gtf format), and extract gene information. Be aware of how the genes are represented.
    """
    f = open(filename, "rt")
    line = f.readline()
    anno_p, anno_n = [], []
    while line:
        if line[0] != '#':
            ls = line.strip().split('\t')
            gtype = ls[2]
            start = ls[3]
            end = ls[4]
            strand = ls[6]
            gene = ls[8]
            if keep_feature and gtype not in keep_feature:
                line = f.readline()
                continue
            gene_name = ''
            for items in gene.strip().split(';'):
                item = items.strip().split('=')
                if item[0] == name:
                    gene_name = item[1]
                    break
            if strand == '+':
                anno_p.append([int(start), int(end), strand, gene_name, gtype])
            else:
                anno_n.append([int(start), int(end), strand, gene_name, gtype])
        line = f.readline()
    f.close()
    return anno_p, anno_n


def iou_search(start, end, anno, cut1=0.5, cut2=0.5, srna='ncRNA'):
    """
    Find genes according to the positions and strands. Be aware of the annotations of the sRNAs.
    """
    iou_max = 0
    pos = -1
    name, gtype = '', ''
    idx = 0
    for i in range(len(anno)):
        iou = (min(end, anno[i][1])-max(start, anno[i][0])) / (max(end, anno[i][1])-min(start, anno[i][0]))
        if iou > 0 and anno[i][4] == 'CDS' and (min(end, anno[i][1])-max(start, anno[i][0]))/(end-start)>=cut1:
            name, gtype =  anno[i][3], anno[i][4]
            idx = 1
        if iou > 0 and anno[i][4] == srna and (min(end, anno[i][1])-max(start, anno[i][0]))/(end-start)>=cut2 and idx != 1:
            name, gtype = anno[i][3], anno[i][4]
            idx = 2
        if iou > iou_max and idx != 1:
            iou_max = iou
            pos = i
    
    if idx != 1 and idx != 2:
        if iou_max > 0:
            name, gtype = anno[pos][3], anno[pos][4]
        else:
            name, gtype = "IGR", "IGR"
    return name, gtype


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Please specify the input files and corresponding features \
    for the accurate annotation step.")

    parser.add_argument('-g', type=str, help='The path of reference genome.')  # reference genome
    parser.add_argument('-t', type=str, help='The path of the annotation file.')  # gff annotation
    parser.add_argument('-i', type=str, help='The path of the interacting regions file.')  # interaction
    parser.add_argument('-c', type=str, help='The path of the chimeric reads file.')  # chimera
    parser.add_argument('-r1', type=str, help='The path of the read_1 fastq file.')
    parser.add_argument('-r2', type=str, help='The path of the read_2 fastq file.')
    parser.add_argument('-b1', type=str, help='The path of the alignment file of read_1.')
    parser.add_argument('-b2', type=str, help='The path of the alignment file of read_2.')
    parser.add_argument('-m', type=str, default='./meme', help='The path of the MEME algorithm implement.')
    parser.add_argument('-d', type=str, default='./motifs', help='The directory for output motifs.')
    parser.add_argument('--srna', type=str, default='ncRNA', help='The name of the sRNA in the annotation file.')
    parser.add_argument('--name', type=str, default='gene', help='The interested feature in the annotation file.')

    args = parser.parse_args()

    # read bam file and store reads information
    reads_dict = defaultdict(set)

    bam1 = pysam.AlignmentFile(args.b1, "rb")
    for read in bam1.fetch():
        reads_dict[read.query_name].add((read.get_blocks()[0][0], read.get_blocks()[-1][1]))

    bam2 = pysam.AlignmentFile(args.b2, "rb")
    for read in bam2.fetch():
        reads_dict[read.query_name].add((read.get_blocks()[0][0], read.get_blocks()[-1][1]))
    print('load bam file')

    # load information of chimeric reads and interaction regions
    interaction = pd.read_table(args.i, sep='\t')
    chimera = pd.read_table(args.c, sep='\t', header=None)

    for i in range(interaction.shape[0]):
        RNA1_start = interaction.iloc[i, 1]
        RNA1_end = interaction.iloc[i, 2]
        RNA1_strand = interaction.iloc[i, 3]
        RNA2_start = interaction.iloc[i, 5]
        RNA2_end = interaction.iloc[i, 6]
        RNA2_strand = interaction.iloc[i, 7]

        reads = chimera.loc[(chimera[1] <= RNA1_end) &
                            (chimera[1] >= RNA1_start) &
                            (chimera[2] == RNA1_strand) &
                            (chimera[4] <= RNA2_end) &
                            (chimera[4] >= RNA2_start) &
                            (chimera[5] == RNA2_strand)]
        for j in range(reads.shape[0]):
            read_name = reads.iloc[j, 6]
            poses = reads_dict[read_name]
            for pos in poses:
                if 0 < RNA1_start - pos[0] <= 100:
                    RNA1_start = pos[0]
                if 0 < RNA2_start - pos[0] <= 100:
                    RNA2_start = pos[0]
                if 0 < pos[1] - RNA1_end <= 100:
                    RNA1_end = pos[1]
                if 0 < pos[1] - RNA2_end <= 100:
                    RNA2_end = pos[1]
        interaction.iloc[i, 1] = RNA1_start
        interaction.iloc[i, 2] = RNA1_end
        interaction.iloc[i, 5] = RNA2_start
        interaction.iloc[i, 6] = RNA2_end
    print('load interaction information')

    # annotation

    keep_feature = ['CDS', 'rRNA', 'pseudogene', 'ncRNA', 'tRNA', 'tmRNA', '6S_RNA', 'RNase_P_RNA', '4.5S_RNA', '5UTR',
                    '3UTR']
    fa = pysam.Fastafile(args.g)

    anno_p, anno_n = read_gtf(args.t, keep_feature, args.name)
    anno = anno_p + anno_n
    gene_list1, gene_list2 = [], []
    type_list1, type_list2 = [], []
    fa1 = []
    fa2 = []
    for i in range(interaction.shape[0]):
        gene1 = ''
        gene2 = ''
        RNA1_start, RNA1_end, RNA1_strand, RNA2_start, RNA2_end, RNA2_strand = \
            interaction.iloc[i, 1], interaction.iloc[i, 2], interaction.iloc[i, 3], interaction.iloc[i, 5], \
            interaction.iloc[i, 6], interaction.iloc[i, 7]
        if RNA1_strand == '+':
            gene1, type1 = iou_search(RNA1_start, RNA1_end, anno_p, 0.5, 0.5, args.srna)
        else:
            gene1, type1 = iou_search(RNA1_start, RNA1_end, anno_n, 0.5, 0.5, args.srna)
        if RNA2_strand == '+':
            gene2, type2 = iou_search(RNA2_start, RNA2_end, anno_p, 0.5, 0.5, args.srna)
        else:
            gene2, type2 = iou_search(RNA2_start, RNA2_end, anno_n, 0.5, 0.5, args.srna)
        #     gene1, type1 = iou_search(RNA1_start, RNA1_end, RNA1_strand, anno, 0.5, 0.5)
        #     gene2, type2 = iou_search(RNA2_start, RNA2_end, RNA2_strand, anno, 0.5, 0.5)
        gene_list1.append(gene1)
        gene_list2.append(gene2)
        type_list1.append(type1)
        type_list2.append(type2)
        try:
            fa1.append(fa.fetch(interaction.iloc[i, 0], interaction.iloc[i, 1], interaction.iloc[i, 2] + 1))
        except:
            fa1.append("")
        try:
            fa2.append(fa.fetch(interaction.iloc[i, 4], interaction.iloc[i, 5], interaction.iloc[i, 6] + 1))
        except:
            fa2.append("")
    interaction['gene1'] = gene_list1
    interaction['gene2'] = gene_list2
    interaction['type1'] = type_list1
    interaction['type2'] = type_list2
    interaction['fasta1'] = fa1
    interaction['fasta2'] = fa2

    # interaction.to_csv(args.i[:-4] + '_anno.csv', index=0)
    interaction.to_csv(f'{args.i[:-4]}_anno.csv', index=0)
    print('annotation')

    # find motifs
    interaction = interaction[(interaction['type2'] == args.srna) & (interaction['type1'] == 'CDS')]
    sRNA_genes = set(interaction['gene2'])

    os.system('mkdir -p ' + args.d)
    os.system('chmod u+x ' + args.m)
    for gene in sRNA_genes:
        df = interaction[interaction['gene2'] == gene]
        if df.shape[0] < 10:
            continue
        fasta_directory = os.path.join(args.d, gene + '.fa')
        out_fa = open(fasta_directory, "wt")
        for i in range(df.shape[0]):
            if df.iloc[i, 2] - df.iloc[i, 1] + 1 <= 8 or df.iloc[i, -2] == '':
                continue
            out_fa.write(">%d_%s\n" % (i + 1, df.iloc[i, -6]))
            # out_fa.write(fa.fetch("NC_000913.3", df.iloc[i, 1] - 1, df.iloc[i, 2]) + '\n')
            out_fa.write(f'{df.iloc[i, -2]}\n')
        out_fa.close()
        # command0 = args.m + ' ' + fasta_directory + ' -oc ' + os.path.join(args.d, gene) + \
        #           ' -dna -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0'
        command = f'{args.m} {fasta_directory} -oc {os.path.join(args.d, gene)} -dna -nostatus -time 14400 -mod zoops -nmotifs 3 -minw 6 -maxw 12 -objfun classic -revcomp -markov_order 0'
        os.system(command)
    print('motifs done')



