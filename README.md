GenomonPipeline
===============
GenomonPipeline, a cancer genome and RNA sequencing data analysis pipeline, efficiently detects of genomic variants and transcriptomic changes. Additionally, it automatically produces rich analysis reports describing data qualities and summary of detected variants.  Users can run GenomonPipeline with ease of use in the HGC supercomputer and also can run it in other high performance computers.

## Manual
http://genomon-project.github.io/GenomonPages/

## Setup

0. Precondition

Make DRMAA and singularity available beforehand.

1. Install

```
virtualenv -p python3 ~/venv/genomon_snakemake
source ~/venv/genomon_snakemake/bin/activate
git clone -b snakemake https://github.com/Genomon-Project/GenomonPipeline.git GenomonPipeline_snakemake
cd GenomonPipeline_snakemake
python setup.py install
pip install snakemake pyyaml drmaa
```

2. Edit sample.csv

Edit pathes of sequence files.
```
vi ./test/5929_sample.csv
```

3. Pull container images

```
singularity pull docker://genomon/bwa_alignment:0.2.0
```

4. Edit config file

Edit `image` options, to pulled `.simg`.
And edit pathes of reference files.
```
vi ./test/dna_genomon.cfg
```

## How to use

1. Configure

```
genomon_pipeline dna ./test/5929_sample.csv ${output_dir} ./test/dna_genomon.cfg
```

2. `snakemake`
```
cd ${output_dir}
snakemake
```

case, dry-run
```
snakemake -n
```

case, re-run (force all)
```
snakemake --forceall
```

case, with dot
```
snakemake --forceall --dag | dot -Tpng > dag.png
```

![](./dag.png)
