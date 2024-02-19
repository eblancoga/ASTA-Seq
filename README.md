#  README of ASTA-Seq (Allele-Specific Targeted Amplicon Sequencing) v 1.0

E. Blanco and M. Barrero (CRG, 2024)

`enrique.blanco at crg.eu`
`mercedes.barrero at crg.eu`

## Table of Contents
1. What is ASTA-Seq?
2. ASTA-Seq Command Line
3. How to Cite

## 1. What is ASTA-Seq

ASTA-Seq (Allele-Specific Targeted Amplicon Sequencing) is an automatical pipeline implemented in two Perl scripts to
extract and summarize information from FASTQ reads containing information about the allele-specific hydroxymethylation
status of regions surrounding several Differentially (hydroxy)Methylated CpGs in X-reactivating gene promoters by DNA
methylation arrays. Regions to be analyzed were selected to include species-specific SNPs and at least one D(h)MP while
maintaining an amplicon size of 200-500 bp. Paired oxidative Bisulfite (OxBS) and the mock-oxidation reaction (Bisulfite,
BS) treated DNA from day 5 reprogramming samples were used as templates for a two-step PCR amplification protocol before
library preparation for high-throughput sequencing.

Analysis of 5mC and 5hmC percentages was done independently for each CpG contained in the PCR amplicons utilizing ASTA-Seq.
Amplicons which contained the SNP and at least one CpG in the same read were analyzed, with coverages ranging from around
14000 to 600000 reads. First, reads were assigned to their corresponding genes by identification of the primer used for the
PCR 2. Then, a 15-20 nucleotides sequence upstream of the SNP was used to classify the reads in Mus or Cas. Next, a 15-20 nucleotides
sequence upstream of the CpG was used to determine the presence of CG or TG. Percentages of CG were calculated independently
for Mus and Cas, in control or IFN gamma-treated samples, in BS or OxBS samples. 5mC percentages were calculated by analyzing
the percentage of CG in OxBS samples, while 5hmC percentages were calculated as %CG (BS) - %CG(OxBS). In BS samples,
unmethylated Cs are converted to Ts, while methylated Cs (5mC and 5hmC) stay as Cs. In OxBS samples, 5hmC is first oxidized
and then converted into Ts, together with unmethylated Cs, while 5mC stays as C. 

See main reference (Barrero et al. 2024) for further details on the generation of the sequencing experiments.

## 2. ASTA-Seq Command Line

The source code consists of two Perl scripts that must be executed sequentially over each raw data file (FASTQ)
to generate the final list of stats per gene, SNP and CG:

### SETUP

#### List of Files and Folders

The latest full distribution release can be downloaded from GitHub:

    https://github.com/eblancoga/ASTA-Seq

This is the current list of archives and folders:

* README:
General description of the software and basic instructions to run ASTA-Seq in your computer.
* LICENSE:
Open software license (GPL version 2.0/3.0).
* data/:
Configuration file used for the samples from the publication.
* scripts/:
Perl scripts to execute the ASTA-Seq pipeline.

#### Configuration file format:

Users will introduce the information about each gene and its associated SNP and GC sequences line by line.
Alternative sequences for each SNP or each GC must be separated by comma. When multiple SNP and/or CG for the
same gene should be quantified separately, every occurrence will be formatted into a different line (and gene
name identifier) of the configuration file. An example of the format that must be followed is:

     #gene   seqF    seqSNP  seqGC   REF     ALT
     ...
     Eif2s3x_4       TAGGGTAGTTTTAGGAAAGGTTTTT       GGGAGTTGGTCGAAA,GGGAGTTGGTTGAAA TATCGATGTGATATTG,TATTGATGTGATATTG       A       G
     ...

### List of commands

Command | Description
--------|-------------
`ASTA-Seq1.pl` | annotate each read with the gene, SNP and CG, and genome from the allele information.
`ASTA-Seq2.pl` | calculate the relative frequencies of each dinucleotide for each gene,CG and species.


### ASTA-Seq step 1. Brief description of the implementation

Using the configuration file information, the script will identify first the gene to which each read belongs.
The input flanking SNP upstream sequence will be used to find the right form of the SNP and deduce the genome (CAS
or MUS) afterwards. Finally, the CG real variation form will be registered and associated to the current read
in process. One read can be assigned only to one gene variant of the configuration file, but during the process
this search will be attempted against all of the primer identifying the same gene. An example of the calling and
the first lines of the output is:

    % perl scripts/ASTA-Seq1.pl data/config_file.txt samples/sample1.fastq > results/COMPACT_sample1.txt

    % cat results/COMPACT_sample1.txt

    ...
    @M03766:461:000000000-L96J3:1:1101:18550:1903   Mtm1_1  GGGAAAATAGAATTTATTGGTTGGTTAGGT  4       TTATGTTAGGTTAAA 229     243     A       GGTTGGTTAGGTTTT 22      36      CG      C,T     A       CAS
    @M03766:461:000000000-L96J3:1:1101:19062:1927   Dlg3_1  GGTTGTTGGGTTAGTAGTGGTGGAG       4       GTCGTTTCGTTTTTT 166     180     G       AAGAGTAGGTTGATT 33      47      CG      G       A       MUS
    @M03766:461:000000000-L96J3:1:1101:13876:1935   Zpf185_1        TTTGTATTAGGTAGGGTTGTTGAGT       5       GTTAGGGAGGTAGT  167     180     G       GTAGGGTTGTTGAGT 15      29      TG      G       A       MUS
    ...

### ASTA-Seq step 2. Brief description of the implementation

This second script takes the result of the first step of the analyzes and calculates the relative frequency of each
combination of gene, species and GC. An exerpt of a call to this script and this class of output files is:

     % perl scripts/ASTA-Seq2.pl -v results/COMPACT_sample1.txt > results/FINAL_sample1.txt
     %%%% 4608397 reads from results/COMPACT_sample1.txt have been acquired 

     % cat results/FINAL_sample1.txt

     ...
     Ddx3x_1 AG      CAS     5       0.0115204718785281      MUS     2       0.00460818875141126     
     Ddx3x_1 CG      CAS     16      0.0368655100112901      MUS     17      0.0391696043869957      
     Ddx3x_1 GG      CAS     4       0.00921637750282252     MUS     3       0.00691228312711689     
     Ddx3x_1 TA      CAS     10      0.0230409437570563      MUS     14      0.0322573212598788      
     Ddx3x_1 TC      CAS     2       0.00460818875141126     
     Ddx3x_1 TG      CAS     21246   48.9527891062418        MUS     22062   50.8329301168176        
     Ddx3x_1 TN      CAS     1       0.00230409437570563     
     Ddx3x_1 TT      CAS     10      0.0230409437570563      MUS     9       0.0207368493813507      
     ...

Full set of samples (raw data files) and final output results are available as processed data at the
GEO record GSE236247 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE236247).

## 3. How to Cite

### Principal Citation

Please, cite the following reference when using ASTA-Seq:

The Interferon Gamma Pathway Enhances Pluripotency and X-Chromosome Reactivation in iPSC Reprogramming.
M. Barrero, A. Lazarenkov, E. Blanco, L.G. Palma, A.V. Lopez-Rubio, M. Bauer, A. Bigas, L. Di Croce, J.L. Sardina and B. Payer.
(under review, 2024).

