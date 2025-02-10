## Dependencies

- data.table: 1.16.4
- ggplot2: 3.5.1
- readxl: 1.4.3
- lme4: 1.1-36
- lmerTest: 3.1-3
- emmeans: 1.10.7
- Seurat: 5.2.1

## Data

Study design and dataset can be found in `data.xlsx` file in the repo. This file is assumed to be under the same folder with all the R scripts.

A public available 10x single-cell PBMC dataset is deposited in the `filtered_gene_bc_matrices` folder.

## How to use

- Each R script can be directly run from command line:

    ``` R
    Rscript [script.R]
    ```

- It can also be run in RStudio.

## Note

- If running in Rstudio, it is recommended to clear memory first.
- I wrote comments to describe what each step does, including my thought process.
- I included a `expect_output` folder in which result tables and plots can be found.

## Q1

`q1.R` script generates:

- `q1.res.tsv`: with sample mean and standard deviation (sd) by trt_group and sample_time for each marker.
- `q1.pdf`: plot sample mean plus/minus 1 unit of sd by trt_group and marker.

## Q2

`q2.R` script generates:

- `q2.res.tsv`: result table with test result.

## Q3

`q3.R` script generates:

- `q3.res.tsv`: result table with test result.
- `q3_plots`: folder where plots are dumped.

## Q4

`q4.R` script generates:

- `q4.rds`: RDS file with `pbmc` Seurat object after filtration based on QC criteria. Should be ready to be loaded using `readRDS()` function.
- `sc_plots`: folder with plots to help guide the decision of QC criteria.
