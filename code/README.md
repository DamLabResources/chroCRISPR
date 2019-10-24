Code
-----

`analysis.ipynb` - Take processed datasets for full analysis and figures presented in manuscript

`rnaseq.ipynb` - Demonstration of gene region file from TSS (transcription start site) to pre-defined promoter region (Default: [-1000, +200])

`GUIDEseq_vs_hotspot2` - Demonstration of DNase-Seq analysis using peak-call tools `hotspot2`. Quantification of relationship between DNase-Seq hotspots and CRISPR-induced cleavage sites was done by the distance of closest neighboring DNase-Seq hotspot to CRISPR-induced cleavage sites. 

`GeneDNaseCoverage.ipynb` - Run `production.py` in python notebook instead of terminal command

`preprocess.py` - A Ruffus supervised pipeline for DNase-Seq mapping

`production.py` - Add DNase-Seq RPM column to CRISPR-induced cleavage sites/Gene promoter

`tasks.py` - functions of running command with arguments including `trim_galore`, `bwa`, `sambamba`, `bedtools`.

---------------------