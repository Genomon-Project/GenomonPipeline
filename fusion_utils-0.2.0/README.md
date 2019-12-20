# fusion_utils

## Introduction

`fusion_utils` is a software for processing the result of [`fusionfusion`](https://github.com/Genomon-Project/fusionfusion).

## Dependency

### Python

Python (>= 2.7), `pysam (>= 0.8.1)`,[`annot_utils`](https://github.com/friend1ws/annot_utils) packages.

### Software

[bedtools](http://bedtools.readthedocs.io/en/latest/)

## Install 
```
git clone  https://github.com/friend1ws/fusion_util.git
cd fusion_utils
python setup.py build install
```


## Commands

### comp
Compare results of fusion calls (by fusionfusion, 
[STAR-fusion](https://github.com/STAR-Fusion/STAR-Fusion), 
[Tophat-fusion](http://ccb.jhu.edu/software/tophat/index.shtml), 
MapSplice2, 
[Genomon-Fusion](http://genomon.hgc.jp/rna/)) or 
SV calls (by GenomonSV).

```
fusion_utils comp [-h] [--margin MARGIN]
                       [--sv_margin_major SV_MARGIN_MAJOR]
                       [--sv_margin_minor SV_MARGIN_MINOR]
                       fusion1.txt
                       {fusionfusion,fusionfusion_part,genomonSV,star_fusion,genomon_fusion,mapsplice2,tophat_fusion}
                       fusion2.txt
                       {fusionfusion,fusionfusion_part,genomonSV,star_fusion,genomon_fusion,mapsplice2,tophat_fusion}
                       output.txt
```

### rmdup

Remove putative duplicates from results of fusion calls.

```
fusion_utils rmdup [-h]
                   [--type {fusionfusion,fusionfusion_part,genomonSV,star_fusion,genomon_fusion,mapsplice2,tophat_fusion}]
                   fusion.txt output.txt
```

### filt
Filter out unreliable candidates the results of fusion calls.

```
fusion_utils filt [-h]
                   [--type {fusionfusion,fusionfusion_part,star_fusion,genomon_fusion,mapsplice2,tophat_fusion}]
                   [--grc] [--genome_id {hg19,hg38,mm10}]
                   [--thres threshould] [--filter_same_gene]
                   [--filter_unspliced]
                   fusion.txt output.txt
```
 
