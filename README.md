# Predicting splice-altering variants

## Usage

### Example

```
  nextflow run main.nf --i input.vcf --indels spliceai_indels.vcf.gz --snvs spliceai_snvs.vcf.gz --fa fasta_ref.fa --gtf gtf_ref.gtf --o /results
```
### Execution parameters

| Parameter   | Type       | Description                           |
| :---------- | :--------- | :---------------------------------- |
| `--i` | `file path/expression` | `Input VCF file path or expression` |
| `--indels` | `file path` | `Precomputed SpliceAI indels scores` |
| `--snvs` | `file path` | `Precomputed SpliceAI snvs scores` |
| `--fa` | `file path` | `FASTA reference file` |
| `--o` | `folder path` | `Output folder path` |
| `--t` | `integer` | `CPU threads to use in each SpliceAI task` |

### Input expressions

Expressions can be utilized as described in the [Nextflow documentation](https://www.nextflow.io/docs/latest/process.html#multiple-input-files) to input multiple files at once. Expressions must be enclosed in quotes.

## Results

Predictions are annotated on the the VCF file info field with the general format:

```
SpliceAI=X|X|X|X;SQUIRLS_SCORE=X
```
