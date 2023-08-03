# RIL-seq-interacting-regions
This is a simple tool to calculate the accurate interacting regions for bacteria genomes. To use this tool, you should first install [RIL-seq](https://github.com/asafpr/RILseq) package. You can either include this script into the whole RIL-seq pipeline, or just take the original output file of RIL-seq as the input of this script.

**Usage:**

```
usage: annotation_motif.py [-h] [-g G] [-t T] [-i I] [-c C] [-r1 R1] [-r2 R2] [-b1 B1] [-b2 B2] [-m M] [-d D] [--srna SRNA] [--name NAME]

Please specify the input files and corresponding features for the accurate annotation step.

optional arguments:
  -h, --help   show this help message and exit
  -g G         The path of reference genome.
  -t T         The path of the annotation file.
  -i I         The path of the interacting regions file.
  -c C         The path of the chimeric reads file.
  -r1 R1       The path of the read_1 fastq file.
  -r2 R2       The path of the read_2 fastq file.
  -b1 B1       The path of the alignment file of read_1.
  -b2 B2       The path of the alignment file of read_2.
  -m M         The path of the MEME algorithm implement.
  -d D         The directory for output motifs.
  --srna SRNA  The name of the sRNA in the annotation file.
  --name NAME  The interested feature in the annotation file.

```

We provide an example for the usage of this tool, please see `run_rilseq.sh`. In this example, we will use public data [SRR9096026/SRR9096032](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131520) for analysis. The result file can be found in `result/hfq_lb.inter.interaction_anno.csv` with annotations of all the regions. By contract, the `result/hfq_lb.inter.interaction.txt` is the original output file from the RIL-seq pipeline.

We also apply [MEME](https://meme-suite.org) motif finding algorithms to the binding targets of different sRNAs. The motifs can be found in the `motif` directory. For example, the interacting regions of `ArcZ` gene has the following motif

![logo1](https://github.com/peipp410/RIL-seq-interacting-regions/blob/master/motifs/arcZ/logo_rc1.png)

In conclusion, the original RIL-seq pipeline cannot find accurate interacting regions due to its low reads utilization rate, and we improve the pipeline by mapping those unused sequencing reads to the reference genome and modify the interacting regions according to the mapped positions of these reads.
