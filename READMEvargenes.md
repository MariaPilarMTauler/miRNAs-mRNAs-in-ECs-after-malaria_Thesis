## RNA-Seq Analysis of Var Genes in Plasmodium falciparum HB3 ##

This pipeline was used to analyse the expression of *var* genes based on exon 1 regions in cultured *Plasmodium falciparum* HB3 parasites.

### Workflow Summary ##

1. **Input:**
   - Raw FASTQ files (paired-end, Illumina)
   - Reference genome: `PlasmoDB-61_PfalciparumHB3_Genome.fasta`
   - Gene annotation: `PlasmoDB-61_PfalciparumHB3.gff`

2. **Pipeline Execution:**
   - Run `run_rnaseq_pipeline.sh` to perform:
     - Quality control (`FastQC`)
     - Alignment (`STAR`)
     - Quantification (`featureCounts`)
     - Reporting (`MultiQC`)
     - fasqfiles names need to be updated

3. **Output:**
   - `output/aligned/`: Aligned BAM files
   - `output/qc/`: Quality control reports
   - `output/counts/`: Gene-level raw counts

4. **Post-processing:**
   - Counts converted to TPM in R via `feature counts to TPM.Rproj`
   - Additional filtering to isolate *var* gene exon 1 regions `export_var_exon1_TPM.R`

### Requirements

- STAR
- FastQC
- Subread (featureCounts)
- R + required packages (e.g., `tximport`, `edgeR`, `tidyverse`)

---

> TPM normalisation and exon-1 selection were critical for analysing differential expression of *var* genes between different parasite samples.
