chromCRISPR
===============================

This repository demonstrates a series of analyses to assess the impact of DNA accessibility on CRISPR-Cas9 cleavage efficiency. This analysis takes multiple open access genome-scale datasets including GUIDE-Seq, CIRCLE-Seq, DNase-Seq and RNA-Seq to systematically characterize crucial determinants for CRISPR-induced editing efficiency. 

**Highlighting results:**
* The condensed chromatin conformation has the potential to abrogate the correlation between gRNA:target similarity and CRISPR-induced cleavage frequency.
* CRISPR-induced sequence editing is possible even in regions where the vast majority of endogenous genes are silent.


Quick review
--------
Full analysis and figure generation could be found in [this python notebook](code/analysis.ipynb)

**List of files required to run full analysis:**
* `processed_data/20181216_GUIDE_sup_data_RPM_clean.csv` 
* `processed_data/20181216_CIRCLE_sup_data_RPM_clean.csv`
* `raw_data/transcriptome/HK_genes.txt`
* `raw_data/transcriptome/HEK/paired/SRR3997505/abundance.tsv`
* `raw_data/transcriptome/U2OS/ERR191523_trimmed/ERR191523/abundance.tsv`
* `raw_data/transcriptome/U2OS_990_TSS_1000up_200down_DNaseSeq.csv`
* `raw_data/transcriptome/HEK_TSS_1000up_200down_DNaseSeq.csv`
-------------------
-------------------

Detailed pipeline
------------
The pipelines that generates all required files to run full analysis descrbied above take publically available dataset from multiple resources:

### Preprocessing of DNase-Seq datasets:
Prepare the sorted bam file for the calculation of Read count Per Million mappable reads per basepair (RPM).

| Dataset | Assay | Link |
| ---     | ---   | ---  |
| HEK293T DNA accessibility | DNase-Seq | [ENCFF774HUB.bam](https://www.encodeproject.org/files/ENCFF774HUB/@@download/ENCFF774HUB.bam) |
| U2OS DNA accessibility | DNase-Seq | [SRR4413990.fastq](https://www.ncbi.nlm.nih.gov/sra/?term=SRR4413990) |

Preprocessing of DNase-Seq on HEK293T:

`python code/preprocessing.py -i raw_data/HEK/ENCFF774HUB.bam -p 'bam' -r raw_data/HG19.fasta`

Preprocessing of DNase-Seq on U2OS:

`python code/preprocessing.py -i raw_data/U2OS/SRR4413990.fastq.gz -p 'fastq' -r raw_data/HG19.fasta`

--------------
*Note:* 

*HG19 reference genome could be downloaded from:* 
[GRCh37/hg19](ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/ARCHIVE/BUILD.37.3/Assembled_chromosomes/chr_accessions_HuRef)

*and indxed by:*

`bwa index -a bwtsw raw_data/HG19.fasta`

--------------

### Calculate DNase-Seq RPM to CRISPR-induced cleavage sites/Gene promoter:
| Dataset | Assay/Platform | Link |
| ---     | ---   | ---  |
| GUIDE-Seq identified | GUIDE-Seq | [GUIDEseq_allgRNAs_identified](raw_data/GUIDEseq_allgRNAs_identified.csv) or [Supplementary Table 2](https://media-nature-com.ezproxy2.library.drexel.edu/original/nature-assets/nbt/journal/v33/n2/extref/nbt.3117-S2.xlsx) |
| CIRCLE-Seq identified | CIRCLE-Seq | [CIRCLEseq_allgRNAs_identified](raw_data/CIRCLEseq_allgRNAs_identified_matched.csv) or [Supplementary Table 2](https://media-nature-com.ezproxy2.library.drexel.edu/original/nature-assets/nmeth/journal/v14/n6/extref/nmeth.4278-S2.xlsx) |
| HG19 gene coordinates | NCBI RefSeq | [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=698225799_VBVYkJAZFNjxKaJnafKIkQY9ZcPB&boolshad.hgta_printCustomTrackHeaders=0&hgta_ctName=tb_ncbiRefSeq&hgta_ctDesc=table+browser+query+on+ncbiRefSeq&hgta_ctVis=pack&hgta_ctUrl=&fbQual=whole&fbUpBases=200&fbExonBases=0&fbIntronBases=0&fbDownBases=200&hgta_doGetBed=get+BED) |
| HEK293T transcriptome | NextSeq 500 | [SRR3997505](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3997505) | 
| U2OS transcriptome | HighSeq 2000 | [ERR191523](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR191523) |

_For pre-defined promoter regions, use this file:_
[hg19_allTSS_1000up_200down.bed.gz](raw_data/transcriptome/hg19_allTSS_1000up_200down.bed.gz)

_For more detail of gene region file, please see:_ [rnaseq.ipynb](code/rnaseq.ipynb)

Add DNase-Seq RPM column to CRISPR-induced cleavage sites/Gene promoter:

e.g.

DNase-Seq RPM on HEK293T gene promoter:

`python code/production.py -L raw_data/transcriptome/hg19_allTSS_1000up_200down.bed.gz -c "HEK293T" -b raw_data/HEK/HEK.se50.DNaseSeq.sorted.bam -o processed_data/HEK_TSS_1000up_200down_DNaseSeq.csv`

DNase-Seq RPM on U2OS gene promoter:

`python code/production.py -L raw_data/transcriptome/hg19_allTSS_1000up_200down.bed.gz -c "U2OS" -b raw_data/U2OS/SRR4413990_trimmed.sorted.bam -o processed_data/U2OS_TSS_1000up_200down_DNaseSeq.csv`

DNase-Seq RPM on HEK293T GUIDE-Seq identified cleavage sites:

`python code/production.py -L raw_data/GUIDEseq_allgRNAs_identified.csv -c "HEK293T" -b raw_data/HEK/HEK.se50.DNaseSeq.sorted.bam -w 100 -o processed_data/HEK_GUIDESeq_DNaseSeq.csv`

DNase-Seq RPM on U2OS GUIDE-Seq identified cleavage sites:

`python code/production.py -L raw_data/GUIDEseq_allgRNAs_identified.csv -c "U2OS" -b raw_data/U2OS/SRR4413990_trimmed.sorted.bam -w 100 -o processed_data/U2OS_GUIDESeq_DNaseSeq.csv`

DNase-Seq RPM on HEK293T CIRCLE-Seq identified cleavage sites:

`python code/production.py -L raw_data/CIRCLEseq_allgRNAs_identified_matched.csv -c "HEK293T" -b raw_data/HEK/HEK.se50.DNaseSeq.sorted.bam -w 100 -o processed_data/HEK_CIRCLESeq_DNaseSeq.csv`

DNase-Seq RPM on U2OS CIRCLE-Seq identified cleavage sites:

`python code/production.py -L raw_data/CIRCLEseq_allgRNAs_identified_matched.csv -c "U2OS" -b raw_data/U2OS/SRR4413990_trimmed.sorted.bam -w 100 -o processed_data/U2OS_CIRCLESeq_DNaseSeq.csv`

----------------

### Ready for figure generation:

Open [analysis.ipynb](code/analysis.ipynb) and displace the correct file names in corresponding input DataFrames.


