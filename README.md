# Predicting splice-altering variants

## Usage

### Example

```
  nextflow run main.nf --i input.vcf --indels spliceai_indels.vcf.gz --snvs spliceai_snvs.vcf.gz --fa fasta_ref.fa --sdb squirls_db --t n --o /results
```
### Execution parameters

| Parameter   | Type       | Description                           |
| :---------- | :--------- | :---------------------------------- |
| `--i` | `file path/expression` | `Input VCF file path or expression` |
| `--indels` | `file path` | `Precomputed SpliceAI indels scores` |
| `--snvs` | `file path` | `Precomputed SpliceAI snvs scores` |
| `--sdb` | `folder path` | `SQUIRLS database folder` |
| `--fa` | `file path` | `FASTA reference file` |
| `--t` | `integer` | `CPU threads to use in each SpliceAI task` |
| `--o` | `folder path` | `Output folder path` |


### Input expressions

Expressions can be utilized as described in the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#multiple-input-files) to input multiple files at once. Expressions must be enclosed in quotes.

## Results

Predictions are annotated on the the VCF file info field with the general format:

```
SQUIRLS_SCORE=X;SpliceAI=X|X|X|X
```

## References

Jaganathan K, Kyriazopoulou Panagiotopoulou S, McRae JF, et al. Predicting Splicing from Primary Sequence with Deep Learning. Cell. 2019;176(3):535-548.e24. doi:10.1016/j.cell.2018.12.015

Danis D, Jacobsen JOB, Carmody LC, et al. Interpretable prioritization of splice variants in diagnostic next-generation sequencing [published correction appears in Am J Hum Genet. 2021 Nov 4;108(11):2205]. Am J Hum Genet. 2021;108(9):1564-1577. doi:10.1016/j.ajhg.2021.06.014