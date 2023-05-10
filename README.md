# Predicting splice-altering variants

## Usage

### Example

```
  nextflow run main.nf --i sample.vcf --pcv 1 --indels spliceai_scores.raw.indel.hg38.vcf.gz --snvs spliceai_scores.raw.snv.hg38.vcf.gz --sdb 2203_hg38 --ref grch38 --fa hg38.fa.gz --t 10 --o results/ -resume
```
### Execution parameters

| Parameter   | Type       | Description                           |
| :---------- | :--------- | :---------------------------------- |
| `--i` | `file path/expression` | `Input VCF file path or expression` |
| `--pcv` | `integer` | `1 if using precalculated scores, 0 if not` |
| `--indels` | `file path` | `Precomputed SpliceAI indels scores` |
| `--snvs` | `file path` | `Precomputed SpliceAI snvs scores` |
| `--sdb` | `folder path` | `SQUIRLS database folder` |
| `--fa` | `file path` | `FASTA reference file` |
| `--ref` | `string` | `Gene annotation file (either 'grch37', 'grch38' or a custom file)` |
| `--t` | `integer` | `CPU threads to use in multithreading-compatible tasks` |
| `--o` | `folder path` | `Output folder path` |


### SQUIRLS database folder

SQUIRLS's database can be downloaded from it's [documentation page](https://squirls.readthedocs.io/en/master/setup.html#squirls-downloadable-resources). Make sure to download the adequate files according to the genome build being used.

### Illumina's SpliceAI precalculated scores

Illumina made available annotations for all possible substitutions, 1 base insertions, and 1-4 base deletions within genes for download at their [online platform](https://basespace.illumina.com/s/otSPW8hnhaZR). These annotations are free for academic and not-for-profit use; other use requires a commercial license from Illumina, Inc. To use them in our pipeline, you should execute it with `pcv = 1` and use the parameters `indels` and `snvs` to indicate the score file paths.

### Input expressions

To input multiple files at once, expressions can be utilized as described in the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#multiple-input-files). Expressions must be enclosed in quotes while individual file paths do not.

### Custom gene annotation files

It's possible to create custom gene annotation files using the files [here](https://github.com/Illumina/SpliceAI/tree/master/spliceai/annotations) as a template.

### CPU threads

Some tasks in the pipeline allow for the usage of multithreading, and the amount of CPU threads used for those will be determined by the `--t` parameter

## Results

Predictions are annotated on the the VCF file info field with the general format:

```
SQUIRLS_SCORE=SCORE;SpliceAI=ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL;SplicePrediction=PREDICTION
```

### SQUIRLS

SQUIRLS provides a single score that is the maximum predicted splicing pathogenicity score among the ones calculated for each variant. More information can be found at the [SQUIRLS documentation](https://squirls.readthedocs.io/en/master/interpretation.html)

### SpliceAI

SpliceAI provides a series of scores and information separated by `|`, in the order presented in the table below. More information can be found at the [SpliceAI Github repository](https://github.com/Illumina/SpliceAI)

| Field   | Description       |
| :---------- | :--------- |
| `ALLELE` | `Alternate allele` |
| `SYMBOL` |  `Gene symbol` |
| `DS_AG` | `Delta score (acceptor gain)` |
| `DS_AL` | `Delta score (acceptor loss)` |
| `DS_DG` | `Delta score (donor gain)` |
| `DS_DL` | `Delta score (donor loss)` |
| `DP_AG` | `Delta position (acceptor gain)` |
| `DP_AL` | `Delta position (acceptor loss)` |
| `DP_DG` | `Delta position (donor gain)` |
| `DP_DL` | `Delta position (donor loss)` |

### Prediction

The `SplicePrediction` is the final annotation of the pipeline for each variant according to the cutoffs and annotation mode utilized. It indicates `P` for a pathogenic prediction and `N` for a non-pathogenic prediction.

## References

Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015

Danis D, Jacobsen JOB, Carmody LC, et al. Interpretable prioritization of splice variants in diagnostic next-generation sequencing [published correction appears in Am J Hum Genet. 2021 Nov 4;108(11):2205]. Am J Hum Genet. 2021;108(9):1564-1577. doi:10.1016/j.ajhg.2021.06.014