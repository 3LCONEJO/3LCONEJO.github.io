---
title: "Day 2: VEP and more bcftools analyses"
date: "2025-02-10"
author: "Jorge Alfredo Suazo Victoria"
output: html_document
---

## 1. Thinking about our experiment

{{< code language="bash" >}}

bcftools isec -C SRR445716.vcf.gz SRR445715.vcf.gz \>\
present_in_IMW004_absent_in_CEN.PK113-7D.txt

{{< /code >}}

```         
chrI    244 C   CT  10
chrI    675 A   G   10
chrI    1152    T   G   10
chrI    1397    A   G   10
chrI    1428    T   C   10
chrI    1757    G   T   10
chrI    2002    G   T   10
chrI    2029    T   C   10
chrI    2406    A   C   10
chrI    12227   C   T   10
```

```         
About:   Create intersections, unions and complements of VCF files.
Usage:   bcftools isec [options] <A.vcf.gz> <B.vcf.gz> [...]

Options:
    -c, --collapse STRING          Treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]
    -C, --complement               Output positions present only in the first file but missing in the others
    -e, --exclude EXPR             Exclude sites for which the expression is true
    -f, --apply-filters LIST       Require at least one of the listed FILTER strings (e.g. "PASS,.")
    -i, --include EXPR             Include only sites for which the expression is true
        --no-version               Do not append version and command line to the header
    -n, --nfiles [+-=~]INT         Output positions present in this many (=), this many or more (+), this many or fewer (-), the exact (~) files
    -o, --output FILE              Write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -p, --prefix DIR               If given, subset each of the input files accordingly, see also -w
    -r, --regions REGION           Restrict to comma-separated list of regions
    -R, --regions-file FILE        Restrict to regions listed in a file
        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -t, --targets REGION           Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT              Use multithreading with <int> worker threads [0]
    -w, --write LIST               List of files to write with -p given as 1-based indexes. By default, all files are written

Examples:
   # Create intersection and complements of two sets saving the output in dir/*
   bcftools isec A.vcf.gz B.vcf.gz -p dir

   # Filter sites in A and B (but not in C) and create intersection
   bcftools isec -e'MAF<0.01' -i'dbSNP=1' -e - A.vcf.gz B.vcf.gz C.vcf.gz -p dir

   # Extract and write records from A shared by both A and B using exact allele match
   bcftools isec A.vcf.gz B.vcf.gz -p dir -n =2 -w 1

   # Extract and write records from C found in A and C but not in B
   bcftools isec A.vcf.gz B.vcf.gz C.vcf.gz -p dir -n~101 -w 3

   # Extract records private to A or B comparing by position only
   bcftools isec A.vcf.gz B.vcf.gz -p dir -n -1 -c all
```

## Question 17: Can you think of a way to obtain a list of candidates that may underlie the ability of these strains to grow on lactate? Hint: You can assume that variants shared by both IMW004 and IMW005 are likely to have arisen before the start of the experiment (i.e., from the unsequenced initial jen1 delta strain), and therefore are not biologically interesting. How many variants (unfiltered) are in IMW004 that are not shared by any other strain?

### My guess

{{< code language="bash" >}}
bcftools isec -C SRR445716.vcf.gz SRR445715.vcf.gz SRR445717.vcf.gz --output-type v -o IMW004_unique.vcf -w 1

bcftools isec -C SRR445717.vcf.gz SRR445715.vcf.gz SRR445716.vcf.gz --output-type v -o IMW005_unique.vcf -w 1

bgzip IMW004_unique.vcf

bgzip IMW005_unique.vcf

bcftools index IMW004_unique.vcf.gz  

bcftools index IMW005_unique.vcf.gz  

bcftools merge IMW004_unique.vcf.gz IMW005_unique.vcf.gz -o Lac_Uniques

{{< /code >}}

## Question 18: How many variants remain in IMW004 after filtering?

{{< code language="bash" >}}
bcftools filter -i'QUAL>=30 && AD[*:1]>=50 && type="snp"' IMW004_unique.vcf.gz -o IMW004.flt.vcf

bcftools view -H IMW004.flt.vcf | wc -l

{{< /code >}}

`25`

## Question 19: How many variants remain in IMW005 after filtering?

{{< code language="bash" >}}
bcftools filter -i'QUAL>=30 && AD[*:1]>=50 && type="snp"' IMW005_unique.vcf.gz -o IMW005.flt.vcf

bcftools view -H IMW005.flt.vcf | wc -l

{{< /code >}}

`6`

## Question 20: What do all the options that we added to the command mean? Hint: Look at the full options in <http://www.ensembl.org/info/docs/tools/vep/script/vep_options.html>.

-   **`--cache`**:\
    This option tells VEP to use locally cached data for annotation. Caching speeds up the annotation process by avoiding repeated queries to Ensembl's online database.

-   **`--dir_cache /home/drobles/.vep/`**:\
    Specifies the directory where the cached files for VEP are stored. In this case, it is pointing to `/home/drobles/.vep/`.

-   **`-i SRR445716_unique.flt.vcf`**:\
    Specifies the input file for VEP. Here, the input is a VCF file named `SRR445716_unique.flt.vcf`, which contains the variants to be annotated.

-   **`-o SRR445716_unique.flt.vep.vcf`**:\
    Specifies the output file name. VEP will write the annotated variants to `SRR445716_unique.flt.vep.vcf`.

-   **`--vcf`**:\
    This option tells VEP to produce output in VCF format. The annotated variants will be written as an updated VCF file.

-   **`--species "saccharomyces_cerevisiae"`**:\
    Specifies the species to be used for annotation. In this case, the annotation will be done for *Saccharomyces cerevisiae* (yeast).

## Question 21: Look at the output VCF. What happened to the original VCF? Did VEP add an annotation? Which one?

Yes, VEP added an annotation to the original VCF. They added the following header

{{< code language="bash">}}
INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID"

{{< /code >}}

# 3. Making sense of our results

Run VEP on both IMW004 and IMW005 filtered VCF files. Study the output very well. Now create a program in your favourite language that outputs:

- Genes are mutated in any or both of the files

- What mutation is present in what strain

Did you find the original mutations found by the authors in the ADY2 gene?

Question 22: Filter the consequences to only keep those that are either missense, stop gained, frameshift, splice acceptor or splice donor. These are typically the mutations that are predicted to directly affect protein function. How many genes are mutated with any of these consequences in both strains?

{{< code language="r" >}}
library(vcfR)
library(tidyverse)

vcf_file <- "SRR445716_unique.flt.vep.vcf"  # Reemplaza con la ruta de tu archivo [IMW004]
vcf <- read.vcfR(vcf_file)
vcf_data <- as.data.frame(vcf@fix)  # Información básica de las variantes
info_data <- vcfR::extract_info_tidy(vcf) 

vcf_data <- vcf_data[-8]

final_data <- vcf_data %>%
  bind_cols(info_data) 

# Divide las anotaciones de INFO, enfocándote en el campo CSQ (anotaciones de VEP)
csq16_data <- final_data %>%
  separate_rows(CSQ, sep = ",") %>%
  separate(CSQ, into = c("Allele", "Consequence", "Impact", "Gene", "Feature", "Feature_type",
                         "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position",
                         "CDS_position", "Protein_position", "Amino_acids", "Codons",
                         "Existing_variation", "Distance", "STRAND", "SYMBOL", "SYMBOL_SOURCE",
                         "HGNC_ID", paste0("Extra_", 23:50)),  # Agregar más columnas para los valores extra
           sep = "\\|", fill = "right")


vcf_file <- "SRR445717_unique.flt.vep.vcf"  # Reemplaza con la ruta de tu archivo [IMW005]
vcf <- read.vcfR(vcf_file)
vcf_data <- as.data.frame(vcf@fix)  # Información básica de las variantes
info_data <- vcfR::extract_info_tidy(vcf) 

vcf_data <- vcf_data[-8]

final_data <- vcf_data %>%
  bind_cols(info_data) 

# Divide las anotaciones de INFO, enfocándote en el campo CSQ (anotaciones de VEP)
csq17_data <- final_data %>%
  separate_rows(CSQ, sep = ",") %>%
  separate(CSQ, into = c("Allele", "Consequence", "Impact", "Gene", "Feature", "Feature_type",
                         "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "cDNA_position",
                         "CDS_position", "Protein_position", "Amino_acids", "Codons",
                         "Existing_variation", "Distance", "STRAND", "SYMBOL", "SYMBOL_SOURCE",
                         "HGNC_ID", paste0("Extra_", 23:50)),  # Agregar más columnas para los valores extra
           sep = "\\|", fill = "right")


# Encuentra los genes comunes entre csq16_data y csq17_data
common_genes <- intersect(csq16_data$Gene, csq17_data$Gene)

# Mostrar los genes comunes
common_genes

# Filtrar las filas de csq16_data con genes comunes
csq16_common <- csq16_data %>% filter(Gene %in% common_genes)

# Filtrar las filas de csq17_data con genes comunes
csq17_common <- csq17_data %>% filter(Gene %in% common_genes)

# Ver los resultados
csq16_common
csq17_common


csq16_mutations <- csq16_common %>%
  select(POS, REF, ALT, Gene, Allele, Consequence, Impact, Existing_variation, Amino_acids) 


csq17_mutations <- csq17_common %>%
  select(POS, REF, ALT, Gene, Allele, Consequence, Impact, Existing_variation, Amino_acids) 

consequences_of_interest <- c("missense_variant", "stop_gained_variant", "frameshift_variant", "splice acceptor_variant", "splice_donor_variant")

csq16_filtered <- csq16_mutations %>%
  filter(Consequence %in% consequences_of_interest)

csq17_filtered <- csq17_mutations %>%
  filter(Consequence %in% consequences_of_interest)

View(csq16_filtered)

View(csq17_filtered)


{{< /code >}}



## IMW004 mutations


|POS    |REF |ALT |Gene |Allele |Consequence      |Impact   |Existing_variation |Amino_acids |
|:------|:---|:---|:----|:------|:----------------|:--------|:------------------|:-----------|
|132370 |G   |C   |ADY2 |C      |missense_variant |MODERATE |gCt/gGt            |252         |


## IMW005 mutations


|POS    |REF |ALT |Gene |Allele |Consequence      |Impact   |Existing_variation |Amino_acids |
|:------|:---|:---|:----|:------|:----------------|:--------|:------------------|:-----------|
|132470 |G   |C   |ADY2 |C      |missense_variant |MODERATE |Cta/Gta            |219         |
|540930 |C   |G   |     |G      |missense_variant |MODERATE |aCg/aGg            |44          |
