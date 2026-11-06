# TF Screen: Monoclonal EB Support Data

This folder contains supporting data for the TF screen of monoclonal EBs, covering both the arrayed (rep-1) and pooled (rep-2) experiments.

## Files

| File | Description |
|------|-------------|
| `cell_clonotype.rds` | Assignment of cells (`cell_id` column) to individual monoclonal EBs (`clonotype_id` column). |
| `gRNA_matrix.rds` | Read count matrix of gRNAs (rows) by cells (columns). |
| `cell_clonotype_gRNA.rds` | The top-represented gRNA for each monoclonal EB. |

