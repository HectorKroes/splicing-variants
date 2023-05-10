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
| `--t` | `integer` | `CPU threads to use in each SpliceAI task` |
| `--o` | `folder path` | `Output folder path` |


### Precalculated scores

Illumina made available annotations for all possible substitutions, 1 base insertions, and 1-4 base deletions within genes for download at their [online platform](https://basespace.illumina.com/s/otSPW8hnhaZR). These annotations are free for academic and not-for-profit use; other use requires a commercial license from Illumina, Inc. To use them in our pipeline, you should execute it with pcv=1 and use parameters 'indels' and 'snvs' to indicate the score file paths.

### Input expressions

Expressions can be utilized as described in the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#multiple-input-files) to input multiple files at once. Expressions must be enclosed in quotes while individual file paths do not.

### Custom gene annotation files

According to SpliceAI instructions, it's possible to create custom gene annotation files using the files [here](https://github.com/Illumina/SpliceAI/tree/master/spliceai/annotations) as a template.

## Results

Predictions are annotated on the the VCF file info field with the general format:

```
SQUIRLS_SCORE=X;SpliceAI=X|X|X|X
```

## References

Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015

Danis D, Jacobsen JOB, Carmody LC, et al. Interpretable prioritization of splice variants in diagnostic next-generation sequencing [published correction appears in Am J Hum Genet. 2021 Nov 4;108(11):2205]. Am J Hum Genet. 2021;108(9):1564-1577. doi:10.1016/j.ajhg.2021.06.014