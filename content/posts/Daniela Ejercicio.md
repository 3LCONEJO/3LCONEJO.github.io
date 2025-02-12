---
title: "Day One in: Variant calling and Ensembl VEP exercises - LCGEJ"
author: "Jorge Alfredo Suazo Victoria"
date: "2025-02-01"
output: 
  html_document: 
    toc: true
---


``` bash

cd /home/suaria/Documents/variant_calling/data

for sample in SRR445715 SRR445716 SRR445717; do
    samtools stats -r /home/suaria/Documents/variant_calling/data/S288C_ref.fa /home/suaria/Documents/variant_calling/data/${sample}.aligned.sorted.bam > ${sample}.stats
    plot-bamstats -r /home/suaria/Documents/variant_calling/data/other_files/S288C_ref.fa.gc -p ${sample}.graphs/ ${sample}.stats
done
```

## Question 1: What is the percentage of mapped reads in all three files? Check the insert size, GC content, per-base sequence content and quality per cycle graphs. Do they all look reasonable?

The percentage of mapped reads in all three files is:

### **SRR445717**

1.  **Total Reads:** 13,730,526

2.  **Mapped Reads:** 13,230,229 (96.4%)

3.  **Mapped Bases:** 660,554,158 (96.2%)

**Conclusion:**

96.4% of the reads are mapped, which is a very good percentage. This indicates that the majority of the sequences aligned correctly.

**Insert Size:**

The peak is around 250 bp, which is normal for many paired-end sequencing libraries.

**GC Content:**

The graph shows a normal distribution around 50%, which is typical for many species, including humans.

**Per-base Sequence Content:**

There is no sign of deviation in the nucleotide composition, suggesting good sample preparation.

**Quality per Cycle:**

The Phred score is high for the majority of cycles, with a slight drop at the end, which is expected in Illumina sequencing.\

### **SRR445716**

1.  **Total Reads:** 12,870,162

2.  **Mapped Reads:** 12,528,002 (96.6%)

3. **Mapped Bases:** 620,450,305 (96.4%)

**Conclusion:**

It's seem like a good percentage of mapped reads. This indicates that the majority of the sequences aligned correctly.

**Insert Size:**

The peak is around 250 bp, which is normal for many paired-end sequencing libraries.

**GC Content:**

The graph shows a normal distribution around 50%, the peak is at 38.7.

**Per-base Sequence Content:**

There is no sign of deviation in the nucleotide composition, suggesting good sample preparation.

**Quality per Cycle:**

The Phred score is high for the majority of cycles, with a slight drop at the end, which is expected in Illumina sequencing.

### **SRR445715**

1.  **Total Reads:** 17,964,244

2.- **Mapped Reads:** 17,503,811 (97.4%)

3.- **Mapped Bases:** 888,144,619 (96.9%)

**Conclusion:**

97.4% of the reads are mapped, which is a very good percentage. This indicates that the majority of the sequences aligned correctly.

**Insert Size:**

It's seems like it have a slight error at the top of the curve.

**GC Content:**

The graph shows a normal distribution around 50%, the peak is at 40.

**Per-base Sequence Content:**

There is no sign of deviation in the nucleotide composition, suggesting good sample preparation.

**Quality per Cycle:**

There's a slight drop and then a depresion in the middle at graph, suggesting a problem in the quality of the reads.

The error rate in theas sample is higher than the other two samples.

# Generating a pileup


``` bash
samtools mpileup  -f /home/suaria/Documents/variant_calling/data/S288C_ref.fa /home/suaria/Documents/variant_calling/data/SRR445715.aligned.sorted.bam | less -S
```

## Question 2: What is the read depth at position chrI:29519? What is the reference base?
Are they any non-reference bases?

There may be some non-reference bases based on the encoded information, but a more detailed base call breakdown is needed to confirm.

| Chromosome | Position | Reference Base | Depth | Bases |Quality|
|------------|----------|----------------|-------|-------|-------|
| chrI       | 29519    | A              | 56    | ,$,,.,.,..,,,,..,..,..........,.,,.,,........,.........., | BCB=A4BB9>BB@?A>B>B?BBBA@A?CB7?C8AB=@BBBCB=B@@@BBCAC?B00|

## Question 3: What about at position chrI:29522? What is the reference base? Are there any non-reference bases?

There are several lowercase letters in the base call string, specifically: a,

| Chromosome | Position | Reference Base | Depth | Bases |Quality|
|------------|----------|----------------|-------|-------|-------|
| chrI       | 29522    | T              | 46    | aaaaAaAAaAAAAAAAAAAaAaaAaaAAAAAAAAaAAAAAAAAAAa | 8;??>:4BB@BABB;A=BABBCBBB?ABA=CABBBAAABC5CAB00|


# 3. Generating genotype likelihoods and variant calling


``` bash
bcftools mpileup -f /home/suaria/Documents/variant_calling/data/S288C_ref.fa /home/suaria/Documents/variant_calling/data/SRR445715.aligned.sorted.bam | less -S
```

> This is an intermediate output that contains genotype likelihoods (if you don’t remember what this is, go back to your notes on the Bayesian exercises we did!) [`OK I will`]


``` bash
bcftools mpileup -f /home/suaria/Documents/variant_calling/data/S288C_ref.fa /home/suaria/Documents/variant_calling/data/SRR445715.aligned.sorted.bam | bcftools call -m --ploidy 1  | less -S
```

## Question 4: Study the command. Why did we use these settings? If you were performing variant calling in human data, what settings would you use?

THe pipe use the command call and the -m parameter is descripted as "Alternative model for multiallelic and rare-variant calling (conflicts with -c)" and the --ploidy 1 means that the organism is haploid, if we want to use it in humans it would be 2.

## Question 5: What option should we add to only print variant sites?

The option that we should add is -v.


``` bash

for sample in SRR445715 SRR445716 SRR445717; do
bcftools mpileup -a AD -f /home/suaria/Documents/variant_calling/data/S288C_ref.fa /home/suaria/Documents/variant_calling/data/${sample}.aligned.sorted.bam -Ou | bcftools call -mv --ploidy 1 -o ${sample}.vcf
done
```

## Question 6: What is the reference and variant base at position chrIV:122724?
|Reference|Variant|
|---------|-------|
|G        |A      |

## Question 7: What is the total read depth at position chrIV:122724?

 DP=58
   
## Question 8: What is the number of high-quality forward reads supporting the variant call at position chrIV:122724? How many reads support the reference allele?
0
## Question 9: What sort of event is happening at position chrI:29007?
The INDEL classification indicates that the event is a structural variation that involves the insertion of the G base at this position.

## 4.- Variant filtering


``` bash
bcftools query -f'POS = %POS\n' SRR445715.vcf | head
```

```
POS = 83
POS = 136
POS = 137
POS = 139
POS = 262
POS = 286
POS = 305
POS = 457
POS = 476
POS = 485
```


``` bash
bcftools query -f'%POS %REF,%ALT\n' SRR445715.vcf | head
```

```
83 AG,A
136 G,A
137 C,CT
139 TCC,TCCCC
262 A,G
286 A,T
305 C,G
457 CAAA,CAA
476 G,T
485 T,C
```


``` bash
bgzip SRR445715.vcf
bgzip SRR445716.vcf
bgzip SRR445717.vcf
bcftools index SRR445715.vcf.gz
bcftools index SRR445716.vcf.gz
bcftools index SRR445717.vcf.gz
bcftools merge -0 -o combined.vcf SRR445715.vcf.gz SRR445716.vcf.gz SRR445717.vcf.gz
```


``` bash
bcftools query -f'%POS %QUAL [%GT %AD  ] %REF %ALT\n' combined.vcf | head
```
```
83 142.328 1 0,11  1 4,19  1 0,13   AG A
136 148.417 1 0,30  1 1,71  1 0,38   G A
137 110.074 1 10,9  1 23,24  1 6,18   C CT
139 42.4134 1 11,3  1 35,24  1 20,9   T TCC
244 48.0595 0 .  1 9,12  1 4,15   C CT
262 5.85486 1 11,13  0 .  0 .   A G
286 223.417 1 0,43  1 0,88  1 0,42   A T
305 201.416 1 0,50  1 0,69  1 0,37   C G
457 35.4232 1 17,16  1 41,58  0 .   CA C
476 160.421 1 0,32  1 0,44  1 0,24   G T
```

## Question 10: Can you print rows with QUAL bigger than 30 and with at least 50 alternate reads? For this we will need to query the second value of the AD field. Note that the indexes are zero-based; the first AD value is represented as “AD[0]”, therefore the second value must be queried as “AD[1]>=50”. However, you will also need to indicate which sample to look at, to look at any sample you can use the asterisk (e.g. the instruction would look like “AD[*:1]>=50”) Hint: If you get stuck, look at the examples that Petr Danecek (pd3) explained here: https://github.com/samtools/bcftools/issues/757


``` bash
bcftools query -f'%POS %QUAL [%GT %AD  ] %REF %ALT\n' -i'QUAL>30 && AD[*:1]>=50' combined.vcf | head
```

```
36 148.417 1 1,71   G A
286 223.417 1 0,88   A T
305 201.416 1 0,50  1 0,69   C G
457 35.4232 1 41,58   CA C
610 225.417 1 0,58  1 0,106  1 0,77   G A
633 225.168 1 33,102  1 20,69   T C
681 225.417 1 0,69  1 0,51   G A
686 185.809 1 21,66  1 13,51   A G
778 228.323 1 0,68  1 12,63  1 4,53   A G
1008 225.417 1 0,86  1 0,66   A G
```


``` bash
bcftools stats SRR445715.vcf.gz | grep TSTV | cut -f5
```

```
# TSTV, transitions/transversions:
[5]ts/tv
2.90
```



## Question 11: How does the ts/tv change if you apply the filters above? Use the bcftools stats command with the -i option to include calls with QUAL at least 30 and the number of alternate reads at least 50.



``` bash
bcftools stats -i'QUAL>=30 && AD[*:1]>=50' SRR445715.vcf.gz | grep TSTV | cut -f5
```

### TSTV, transitions/transversions:
[5]ts/tv
3.25

## Question 12: What is the ts/tv of removed sites?


``` bash
bcftools stats -e 'QUAL>=30 && AD[*:1]>=50' SRR445715.vcf.gz | grep TSTV | cut -f5
```

### TSTV, transitions/transversions:
[5]ts/tv
2.25

### Confusion

Can you apply the above filters to produce a final callset and apply the -g and -G options to soft filter variants around indels? [Is this a real question?, I thought that the command below was the answer]


``` bash
bcftools filter -sLowQual -m+ -i'QUAL>=30 && AD[*:1]>=50' -g8 -G10 combined.vcf -o combined.flt.vcf
```

## Question 13: Can you see why position chrI:1101 passed these filters? And why did position chrI:1107 fail the filtering step?


```
chrI    1101    .       G       A       152.416 PASS    VDB=0.316037;SGB=-0.69>
chrI    1107    .       T       C       86.1177 LowQual VDB=0.0031237;SGB=-0.6>
```

> INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">

I sincerely don't understand completely what the VDB means, so i asked to my realy close friend.

<a href="#5-multi-sample-variant-calling"><button>Skip This</button></a>

## CHAT GPT

This line is from the **VCF (Variant Call Format)** file header, and it describes the **INFO** field called `VDB`. Here's a breakdown of what each part means:



---

### **1. Field Components:**

- **`ID=VDB`**  
  - This is the **identifier** of the INFO field. In this case, `VDB` stands for **Variant Distance Bias**.

- **`Number=1`**  
  - This indicates that the field contains **one value** per variant.

- **`Type=Float`**  
  - The value for `VDB` is a **floating-point number** (i.e., a decimal number).

- **`Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)"`**  
  - This describes the purpose of the `VDB` field:
    - **Variant Distance Bias (VDB)** measures how evenly the variant-supporting reads are distributed across the sequence.
    - It’s particularly useful for **filtering out artefacts** in **RNA-seq data**, such as false variants near **splice sites**.
    - **“Bigger is better”** means that a **higher VDB** value suggests the variant is **more likely to be real** (less biased), while a lower VDB might indicate an artefact.

- **`Version="3"`**  
  - This indicates the **version** of the software or algorithm that generated this `VDB` value.

---

### **2. What Is Variant Distance Bias (VDB)?**

- **Context:**  
  In RNA-seq data, reads can align imperfectly around **splice sites** (where exons are joined), leading to **artefacts**—false variants caused by misalignments rather than true biological variation.

- **VDB Calculation:**  
  VDB measures how **evenly distributed** the positions of variant-supporting reads are along the alignment:
  - If reads supporting the variant are **clustered** together (e.g., near a splice site), the VDB will be **low**, suggesting a potential artefact.
  - If reads are **spread out** across the region, the VDB will be **high**, indicating the variant is likely genuine.

- **Filtering Use:**  
  You can use VDB to **filter out false positives** in variant calling:
  - **Low VDB**: May indicate artefacts—consider filtering these out.
  - **High VDB**: More reliable variants.

---

### **3. Practical Example:**

Let’s say you’re filtering your VCF file and want to exclude variants with potential splice-site artefacts. You might set a **VDB threshold** to only keep variants with higher values:

```bash
bcftools view -i 'INFO/VDB > 0.05' your_file.vcf > filtered_file.vcf
```

Here, variants with `VDB > 0.05` are retained, assuming values below that might be artefacts.

---

Let me know if you need more details on how to use this in your filtering or interpretation!

So the position chrI:1101 passed the filter because the VDB value was higher (0.316037), suggesting a more evenly distributed set of variant-supporting reads. In contrast, position chrI:1107 failed the filtering step due to its lower VDB value (0.0031237), indicating a potential artefact or bias in the variant-supporting reads.

# 5. Multi-sample variant calling


## Question 14: There are three BAM files in the original directory /mnt/atgc-d2/bioinfoII/drobles/variant_calling/data/ [sure they are]. Can you modify the command from section 3 to use all three BAM files and only write out variant sites in chromosome I? Write the output to a compressed BCF file called multi.bcf and index the file afterwards.

I needed to do this before but, here's the manual

```
Usage: bcftools mpileup [options] in1.bam [in2.bam [...]]

Input options:
  -6, --illumina1.3+      Quality is in the Illumina-1.3+ encoding
  -A, --count-orphans     Do not discard anomalous read pairs
  -b, --bam-list FILE     List of input BAM filenames, one per line
  -B, --no-BAQ            Disable BAQ (per-Base Alignment Quality)
  -C, --adjust-MQ INT     Adjust mapping quality [0]
  -D, --full-BAQ          Apply BAQ everywhere, not just in problematic regions
  -d, --max-depth INT     Max raw per-file depth; avoids excessive memory usage [250]
  -E, --redo-BAQ          Recalculate BAQ on the fly, ignore existing BQs
  -f, --fasta-ref FILE    Faidx indexed reference sequence file
      --no-reference      Do not require fasta reference file
  -G, --read-groups FILE  Select or exclude read groups listed in the file
  -q, --min-MQ INT        Skip alignments with mapQ smaller than INT [0]
  -Q, --min-BQ INT        Skip bases with baseQ/BAQ smaller than INT [1]
      --max-BQ INT        Limit baseQ/BAQ to no more than INT [60]
      --delta-BQ INT      Use neighbour_qual + INT if less than qual [30]
  -r, --regions REG[,...] Comma separated list of regions in which pileup is generated
  -R, --regions-file FILE Restrict to regions listed in a file
      --ignore-RG         Ignore RG tags (one BAM = one sample)
  --ls, --skip-all-set STR|INT  Skip reads with all of the bits set []
  --ns, --skip-any-set STR|INT  Skip reads with any of the bits set [UNMAP,SECONDARY,QCFAIL,DUP]
  --lu, --skip-all-unset STR|INT  Skip reads with all of the bits unset []
  --nu, --skip-any-unset STR|INT  Skip reads with any of the bits unset []
  -s, --samples LIST      Comma separated list of samples to include
  -S, --samples-file FILE File of samples to include
  -t, --targets REG[,...] Similar to -r but streams rather than index-jumps
  -T, --targets-file FILE Similar to -R but streams rather than index-jumps
  -x, --ignore-overlaps   Disable read-pair overlap detection
      --seed INT          Random number seed used for sampling deep regions [0]

Output options:
  -a, --annotate LIST     Optional tags to output; '?' to list available tags []
  -g, --gvcf INT[,...]    Group non-variant sites into gVCF blocks according
                          To minimum per-sample DP
      --no-version        Do not append version and command line to the header
  -o, --output FILE       Write output to FILE [standard output]
  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;
                          'z' compressed VCF; 'v' uncompressed VCF; 0-9 compression level [v]
  -U, --mwu-u             Use older probability scale for Mann-Whitney U test
      --threads INT       Use multithreading with INT worker threads [0]

SNP/INDEL genotype likelihoods options:
  -X, --config STR        Specify platform specific profiles (see below)
  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [20]
  -F, --gap-frac FLOAT    Minimum fraction of gapped reads [0.05]
  -h, --tandem-qual INT   Coefficient for homopolymer errors [500]
  -I, --skip-indels       Do not perform indel calling
  -L, --max-idepth INT    Maximum per-file depth for INDEL calling [250]
  -m, --min-ireads INT    Minimum number gapped reads for indel candidates [2]
  -M, --max-read-len INT  Maximum length of read to pass to BAQ algorithm [500]
  -o, --open-prob INT     Phred-scaled gap open seq error probability [40]
  -p, --per-sample-mF     Apply -m and -F per-sample for increased sensitivity
  -P, --platforms STR     Comma separated list of platforms for indels [all]
  --ar, --ambig-reads STR   What to do with ambiguous indel reads: drop,incAD,incAD0 [drop]
      --indel-bias FLOAT  Raise to favour recall over precision [1.00]
      --indel-size INT    Approximate maximum indel size considered [110]

Configuration profiles activated with -X, --config:
    1.12:        -Q13 -h100 -m1 -F0.002
    illumina:    [ default values ]
    ont:         -B -Q5 --max-BQ 30 -I [also try eg |bcftools call -P0.01]
    pacbio-ccs:  -D -Q5 --max-BQ 50 -F0.1 -o25 -e1 --delta-BQ 10 -M99999

Notes: Assuming diploid individuals.

Example:
   # See also http://samtools.github.io/bcftools/howtos/variant-calling.html
   bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf


```


``` bash
bcftools mpileup -a AD -f /home/suaria/Documents/variant_calling/data/S288C_ref.fa /home/suaria/Documents/variant_calling/data/SRR445715.aligned.sorted.bam /home/suaria/Documents/variant_calling/data/SRR445716.aligned.sorted.bam /home/suaria/Documents/variant_calling/data/SRR445717.aligned.sorted.bam -r chrI -Ou | bcftools call -mv --ploidy 1 -o multi.bcf
```

<a href="/posts/untitled/"><button>Download The Binary</button></a>

## Question 15: Can you apply the same filters as before? How many sites pass the filters? Write the output to a BCF file called multi.filt.bcf and index the file.


<a href="/multi.filt.bcf"><button>Download The Filtered Binary</button></a>

Manual for `filter`
```
About:   Apply fixed-threshold filters.
Usage:   bcftools filter [options] <in.vcf.gz>

Options:
    -e, --exclude EXPR             Exclude sites for which the expression is true (see man page for details)
    -g, --SnpGap INT[:TYPE]        Filter SNPs within <int> base pairs of an indel (the default) or any combination of indel,mnp,bnd,other,overlap
    -G, --IndelGap INT             Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass
    -i, --include EXPR             Include only sites for which the expression is true (see man page for details
        --mask [^]REGION           Soft filter regions, "^" to negate
    -M, --mask-file [^]FILE        Soft filter regions listed in a file, "^" to negate
        --mask-overlap 0|1|2       Mask if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -m, --mode [+x]                "+": do not replace but add to existing FILTER; "x": reset filters at sites which pass
        --no-version               Do not append version and command line to the header
    -o, --output FILE              Write output to a file [standard output]
    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]
    -r, --regions REGION           Restrict to comma-separated list of regions
    -R, --regions-file FILE        Restrict to regions listed in a file
        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]
    -s, --soft-filter STRING       Annotate FILTER column with <string> or unique filter name ("Filter%d") made up by the program ("+")
    -S, --set-GTs .|0              Set genotypes of failed samples to missing (.) or ref (0)
    -t, --targets REGION           Similar to -r but streams rather than index-jumps
    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps
        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]
        --threads INT              Use multithreading with <int> worker threads [0]
```


``` bash
bcftools filter -i'QUAL>=30 && AD[*:1]>=50 && type="snp"' multi.bcf -o multi.filt.bcf

bcftools view -H multi.filt.bcf | wc -l
```

`811` Tara!

Question 16: What is the ts/tv of the raw calls and of the filtered set?


``` bash
bcftools stats multi.filt.bcf | grep TSTV | cut -f5
```

### TSTV, transitions/transversions:
[5]ts/tv
2.56
