# **BETA**
Not in a final state, please do not use unless you are a tester and know how to use it.

# **Context-Aware Transcript Quantification from Long Read Single-Cell and Spatial Transcriptomics data**
This is a pipeline developed for context-aware transcript discovery and quantification from long read single-cell and spatial transcriptomics data. The pipeline consists of 4 steps, barcode/UMI identification and demultiplexing with [flexiplex](https://davidsongroup.github.io/flexiplex/), genome alignment with minimap 2, and transcript discovery and quantification with [Bambu](https://github.com/GoekeLab/bambu/tree/BambuDev). The final output includes novel transcripts found in the sample and transcript level count matrices for each barcode/coordinate.

### **Content** 
- [Installation](#installation)
- [General Usage](#General-Usage)
- [Release History](#Release-History)
- [Citation](#Citation)
- [Contributors](#Contributors)


### **Installation** 
To run this pipeline, you will need [Nextflow](https://www.nextflow.io/docs/latest/install.html) and [Docker](https://docs.docker.com/engine/install/ubuntu/) (or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html_) if you do not have user rights for docker). The latest version is recommended. Note that you mad need to had java to your path if it is not already there.

### **General Usage** 
The pipeline can be run with one line, providing either the path directly to your input fastq/bam file or a samplesheet when running multiple samples. Below are example use cases which you can test using example data that comes with the container.

You may need to specify -r main depending on your local environment. 
Currently investigating issues when running on high performance clusters and pathing to the container.
Note that in the examples that while we do not provide a whitelist paramter, it is generally recommended to do so. See arguments below.

**Running a single sample**
``` 
nextflow run GoekeLab/bambu-SingleCell-Spatial \
  --reads $PWD/examples/reads_chr9_1_1000000.fastq.gz \
  --genome $PWD/examples/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation $PWD/examples/Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \
  --chemistry 10x5v2 \
  --ncore 16 --outdir output \
  -with-singularity lingminhao/bambusc:beta
``` 

**Running multiple samples**
``` 
nextflow run GoekeLab/bambu-SingleCell-Spatial \
  --reads $PWD/examples/samplesheet_basic.csv \   # See the arguments section for format specifications
  --genome $PWD/examples/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation $PWD/examples/Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \
  --chemistry 10x5v2 \ #can be provided in the samplesheet instead
  --ncore 16 --outdir output \
  -with-singularity lingminhao/bambusc:beta
``` 
examples/samplesheet_basic.csv 
| sample            | fastq     |
|:---|:----------|
| replicate1      | reads_chr9_1_1000000.fastq.gz |
| replicate2      | reads_chr9_1_1000000_rep2.fastq.gz |

**Running from a demultiplexed bam (skip demultiplexing and alignment)**
``` 
nextflow run GoekeLab/bambu-SingleCell-Spatial \
  --bams $PWD/examples/demultiplexed.bam \   # See the arguments section for demutiplex format requirements
  --genome $PWD/examples/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation $PWD/examples/Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \ 
  --ncore 16 --outdir output \
  -with-singularity lingminhao/bambusc:beta
``` 

You can run the code above with your own dataset by replacing the arguments accordingly. The pipeline works on both the **long read sequencing** platform Oxford Nanopore Technologies (ONT) and PacBio at the single-cell and spatial level. The description for each argument is shown below: 

### **Required arguments**

The BambuSC pipeline can be started from from raw reads (fastq.gz) or from a demultiplexed .bam file if you have already produced these (from earlier runs of this pipeline or other upstream tools). Therefore either the --reads or --bams argument is mandatory depending on your input files.

**--reads** [PATH] - a path to one read file (.fastq.gz) from a single-cell or spatial experiment run. 

Multiple Samples:

If running multiple samples you can instead path to a csv samplesheet. Samples with the same 'sample' name will be combined together and treated as one sample (for example when combining technical replicates that have the same barcodes). See the following argument descriptions on how the samplesheet can be filled out. 

Running Multiple Similiar Samples 
| sample            | fastq     |
|:---|:----------|
| replicate1      | path/to/1.fastq.gz |
| replicate2      | path/to/2.fastq.gz |
| replicate3      | path/to/3.run1.fastq.gz |
| replicate3      | path/to/3.run2.fastq.gz |

Additionally if your multiple samples require different parameters for demultiplexing and alignment you can provide these in the samplesheet. If any of the chemistry, technology, or whitelist columns are not provided, or if there is a missing entry, the pipeline will use the input from the --chemistry --technology and --whitelist arguments respectively for the samples. 

Multiple Sample with different parameters example: (examples/samplesheet_custom_example.csv)
| sample            | fastq     | chemistry      | technology | whitelist |
|:---|:----------|:------|:------|:------| 
| 3prime_1      | path/to/1.fastq.gz | 10x3v2 | ONT | path/to/whitelist.gz |
| 3prime_2      | path/to/2.fastq.gz | 10x3v2 | ONT | path/to/whitelist.gz |
| 3prime_3_PacBio      | path/to/3.fastq.gz | 10x3v2 | PacBio | path/to/whitelist.gz |
| customSample1      | path/to/4.fastq.gz | -x CTACACGACGCTCTTCCGATCT -b ???????????? -u ???????????? -x TTTTTTTTT | ONT | path/to/whitelist_custom.gz |
| customSample2      | path/to/5.fastq.gz | -x CTACACGACGCTCTTCCGATCT -b ???????????? -u ???????????? -x TTTTTTTTT | ONT | path/to/whitelist_custom.gz |
| customSample3_PacBio      | path/to/6.fastq.gz | -x CTACACGACGCTCTTCCGATCT -b ???????????? -u ???????????? -x TTTTTTTTT | PacBio | path/to/whitelist_custom.gz |

**--bams** [PATH] - If you did demultiplexing and alignment seperately, you can start the pipeline at the bambu step. Requires a path to the alignment file (.bam). The bam file must be demultiplexed in one of the following 3 ways:

1. The read names have the barcode and UMI sequences appended to the start of the name: e.g. `ATCCGTCCAACGGGTA_TTGCTGGCGTGT#46f0ce76-6a12-4a12-a707-2ffef7be7594_+1of1`. The first set of characters (`ATCCGTCCAACGGGTA`) followed by an '_' refer to the cell barcode (CB) and the next set of characters after the underscore  (`TTGCTGGCGTGT`) followed by a '#' refers to the unique molecular identifier (UMI). If there is no UMI in your data this set can be empty and no UMI deduplication will be performed, but ensure the underscore and hash are still present so that the bam file name can be correctly parsed.

2. The barcode is located in the BC tag in the bam file, and the UMI is located in the UG tag

3. The --barcode_map argument is provided in (or as a column in the samplesheet) to a .tsv or .csv with three columns. All reads not present in the table will be discarded from the analysis. The UMI column is optional if you do not need UMI deduplication or do not have UMIs in your data.

Example: NO HEADER IN FILE (see examples/barcode_map_example.csv)

| Read Name            | Barcode     | UMI      | 
|:---|:----------|:------| 
| #46f0ce76-6a12-4a12-a707-2ffef7be7594      | ATCCGTCCAACGGGTA | TTGCTGGCGTGT |

Multiple Samples:

If running multiple samples you can instead path to a csv sample sheet containing the parameters neede for each sample. Samples with the same 'sample' name will be combined together and treated as one sample (for example when combining technical replicates that have the same barcodes). See the following argument descriptions on how the samplesheet can be filled out. Providing a barcode map is optional, and only required if the bam file does not contain the barcode information. Providing a spatial_whitelist is only required if using spatial data, and each sample uses a different whitelist/coordinates, otherwise it will use the whitelist providied by the --whitelist argument.

Multiple Sample Example (examples/samplesheet_bam_example.csv):
| sample            | bam     | barcode_map | spatial_whitelist
|:---|:----------|:------|:------|  
| condition1      | path/to/demultiplexed.bam | | |
| condition2      | path/to/not_demultiplexed.bam | path/to/condition2_barcode_map.csv | |

**--genome** [PATH] - a path to the genome fasta file. in the event you are running the pipeline using a bam file, this should be the same file used for read alignment.

**--annotation** [PATH] - a path to the reference annotation file (.gtf), a txdb object, or annotations object prepared by prepareAnnotations(). When not provided, de novo transcript discovery is performed. (see [Bambu](https://github.com/GoekeLab/bambu) for more information)

**--chemistry** [10x3v2/10x3v3/10x5v2/STRING] - type of 10X Genomics chemistry used in the single-cell or spatial experiment

| Choice of argument            | Description     | 
|:---|:----------------| 
| 10x3v2      | 10x 3' version 2|
| 10x3v3      | 10x 3' version 3|
| 10x5v2      | 10x 5' version 2|

Advanced use: If your sample was not barcoded by one of the above chemistries you can provide a custom search string as detailed in the flexiplex documentation (https://davidsongroup.github.io/flexiplex/). e.g --chemistry '-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTTTTTTT'

### **Output** ###

After the run, all the outputs will be stored in the folder specified by the `--outdir` parameter. The output files are described below and are in 3 categories: General, Single Cell Level, Pseudobulk/EM Level
 
 General
| Output file name                | Description                                                             |
|:----------------------|:------------------------------------------|
| extended_annotations.gtf        | Extended transcript & gene annotations for the genome using long reads data.        |
| txANDGenes.tsv                 | Gene ID associated to each transcript arranged as in the transcript count estimates          |
| genes.tsv                 | Gene ID arranged as in the gene count estimates          |
| extended_annotations_NDR1.gtf                | A gtf file containing all potential transcript models and their NDR score from Bambu (maximum sensitivity).   |

Single-Cell Level
| Output file name                | Description                                                             |
|:----------------------|:------------------------------------------|
| counts_transcript.mtx           | Total read counts estimates for each transcript in each sample (sparse matrix format).        |
| counts_gene.mtx                 | Gene read counts estimates for each transcript in each sample.         |
| incompatibleCounts.mtx                 | The counts that are unable to be assigned to transcripts per gene.         |
| sampleData.tsv                 | Information for the columns for each sparse matrix. Includes if provided the sample names, barcodes, and x and y coordinates (for spatial)         |

Pseudobulk/EM Level
If clusters = "none" then this was also be at the single-cell level. See clusters
| Output file name                | Description                                                             |
|:----------------------|:------------------------------------------|
| counts_transcript.mtx           | Total read counts estimates for each transcript in each sample (sparse matrix format).        |
| CPM_transcript.mtx              | Counts per million (CPM) estimates for each transcript in each sample (sparse matrix format). |
| fullLengthCounts_transcript.mtx | Full length read counts estimates for each transcript in each sample (sparse matrix format).  |
| uniqueCounts_transcript.mtx                | Unique read counts estimates for each transcript in each sample (sparse matrix format).       |
| counts_gene.mtx                 | Gene read counts estimates for each transcript in each sample.         |
| incompatibleCounts.mtx                 | The counts that are unable to be assigned to transcripts per gene.         |
| sampleData.tsv                 | Information for the columns for each sparse matrix. Includes if provided the sample names, barcodes, and x and y coordinates (for spatial)         |
| cellMixs.rds                 | PLACEHOLDER  |
| clusters.rds                 | PLACEHOLDER |


You may then use these count matrices for downstream analysis using tools like [Seurat](https://satijalab.org/seurat/) or [scanpy](https://www.google.com/search?q=scanpy&oq=scanpy&aqs=chrome..69i57.866j0j7&sourceid=chrome&ie=UTF-8). The following short tutorials demonstrate how to analyse the long read [single-cell]() and [spatial]() data using [Seurat](https://satijalab.org/seurat/). These can also be reloaded into R using importBambu("/path/to/outputDir")

### **Optional arguments**

**--whitelist** [PATH] - a path to the `gzipped` cell/spot barcode whitelist or the name of the barcode list ([single-cell](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-) / [spatial](https://kb.10xgenomics.com/hc/en-us/articles/360041426992-Where-can-I-find-the-Space-Ranger-barcode-whitelist-and-their-coordinates-on-the-slide-)) provided by 10X Genomics. It is recommended to provide the correct whitelist as it is used as a reference to filter against low confidence putative barcodes, however you may also run this pipeline without providing it

**--spatial** - a flag to tell the pipeline this is spatial data. If set, the --whitelist needs to also be set and additionally contain the X and Y coordinates for each barcode in two additional columns (tsv with no header). See #Spatial Analysis below for more details

**--fusionMode** - Will run JAFFAL to identify fusion breakpoints and then bambu to identify and quantify fusion transcripts at the single-cell level

**--flexiplex_e** [INT] - the edit distance used by flexiplex for barcode identification (default: 1)

**--technology** [ONT/PacBio]- the long read sequencing platform used to generate the reads [`ONT` or `PacBio`. default:`ONT`]
This changes the stringency of the alignment based on minimap2 recommendations.

**--NDR** [FLOAT] - the NDR threshold used for Bambu transcript discovery. A value between 0 - 1, with lower values being more precise but less sensitive. Typically a value of 0.1 is suitable for most analyses. If you do not want transcript discovery set this value to 0. (default: automatically determined based on data)

**--noEM** - When provided, the pipeline will not perform the final EM quantification step, and the final outputs will be the unique counts and gene counts summarised experiment objects

**--resolution** [FLOAT] - the resolution used for the default Seurat clustering before the EM when --clusters is not provided

**--clusters** [PATH, default: 'auto'] - a path to a tsv or csv (without a header), where the first column contains the barcodes and the second column contains the cluster name. Clustering will be done across all samples with the same sample names. If --lowMemory is provided, clustering is done per input file, they are not combined. Find an example in examples/barcode_clusters_example.csv. Default is 'auto' which will automatically cluster the cells based on gene expression with limited filtering. For advanced use cases where there is sufficient read counts per cell, --clusters can be set to 'none' to have the EM appled to each cell individually.

| barcode            | cluster_id     |
|:---|:----------|:------| 
| GCGCGATAGCTAACAA | cluster1 |
| TAGTGGTTCCTTTCTC | cluster1 |
| TCTGGAAAGGTGACCA | cluster2 |

**--processByBam** [TRUE/FALSE. default: FALSE] - At the cost of runtime, but reducing the max memory the pipeline will use at one time, the pipeline will perform the bambu steps for each input file seperately. *Important* - because read class construction is done seperately using this mode, this will impact the transcripts that are discovered and the quantification results slightly due to small stocastic differences in junction error correction. 

**--processByChromosome** [TRUE/FALSE. default: TRUE]- At the cost of runtime, but reducing the max memory the pipeline will use at one time, the pipeline will perform the bambu steps for each chromosome seperately seperately. *Important* - because read class construction is done seperately using this mode, this will impact the transcripts that are discovered and the quantification results slightly due to small stocastic differences in junction error correction. 

**--ncore** [INT] -  the number of cores to use in the run (default: `12`)

**--cleanReads** [TRUE/FALSE] - if true, only the first supplimentary/primary alignment closest to the barcode is kept for analysis for each read name (default: 'true'). This reduces the liklihood of chimeric reads with undetected barcodes from leading to alignments being incorrectly assigned, however it means only the first alignment of fusion transcripts will be kept.

**--deduplicateUMIs** [TRUE/FALSE] - if true, the longest alignment from the same barcode with the same UMI will be kept, and all others not included in the analysis. Set this to FALSE if UMIs are not provided in the demultiplexing of the bam (default: 'FALSE') 

**--outdir** [PATH] - the path to the output directory

**--bambuPath** [PATH] - For advanced users. The path to a local version of Bambu which will be used instead of the containerized version

### **Multiple Sample Analysis** ###

If you would like to analysis multiple samples or replicates, please note the following differences in the arguments. 
The input to the --reads/--bams arguments needs to be a .csv file with all columns present. See the arguments section to ensure the formatting of the samplesheet.csv is correct and contains all the required information.
The sample column is used to distinguish different replicates and the sample name will be appended to the barcodes from that sample. Therefore if you would like to combine samples that share barcodes you can assign them the same sample name and their counts will be comined. 
The same --genome  and --annotation argument inputs will be applied to all samples. 
An output file is produced for each sample prefixed with the sample name. 

### **Spatial Analysis** ##

To run spatial data using Bambu, the X and Y coordinate mapping for each barcode to a spot needs to be provided as two additional columns in the whitelist (tab seperated). This file can generally be sourced from the provider of the library preperation. 

Example (**Headers are descriptive only, there should be no headers in the input file!**):
| Barcode            | X coordinate | Y coordinate| 
|:---|:---|:---|
| AAACAACGAATAGTTC | 17 | 1 |
| AAACAAGTATCTCCCA  | 103 | 51 |
| AAACAATCTACTAGCA | 44 | 4 |

### **Fusion Transcript Analysis** ##
Two extra parameters are required when running fusion mode
**--fusionMode**
**--jaffal_ref_dir** [PATH] - This needs to be a path to a JAFFA compatible reference. See JAFFAL documentation (https://github.com/Oshlack/JAFFA/wiki/FAQandTroubleshooting#how-can-i-run-jaffa-with-hg19-or-mm10)
When --fusionMode is provided to Bambu-Pipe, three additional steps are included in the pipeline. The first runs JAFFAL identifying the fusion genes breakpoint regions in the sample from the raw fastq files. The next step uses the identified breakpoints to generate artificial fusion scaffolds, placing the fused genes next to each other on their own scaffold. The final step runs Bambu-Clump using the fusion scaffolds and annotations, and will use the clustering provided by --clusters. 

Outputs:

Fusion mode will produce the same outputs alongside the regular steps with the prefix of _fusion and _fusion_EM. The following fusion specific additional outputs are provided.
fusionGene.fasta - These are the fusion scaffolds which the reads are aligned too
fusion.gtf - These are the annotations used by bambu containing the two fused genes, which map to the fusionGene.fasta
jaffa_results.csv - These are the results from JAFFAL describing the fusion breakpoints used to generate the above files.

### **Additional Information**

UMI correction is done at the barcode level. The longest read for each unique barcode-UMI combination is kept for analysis.

To reduce the impact of chimeric reads, only the first alignment after the barcode will be used for analysis.

### **Release History** 

Beta Release: 2023-May-03

### **Citation**
Bambu
Chen, Y., Sim, A., Wan, Y.K. et al. Context-aware transcript quantification from long-read RNA-seq data with Bambu. Nat Methods (2023). https://doi.org/10.1038/s41592-023-01908-w

Flexiplex
Oliver Cheng, Min Hao Ling, Changqing Wang, Shuyi Wu, Matthew E Ritchie, Jonathan Göke, Noorul Amin, Nadia M Davidson, Flexiplex: a versatile demultiplexer and search tool for omics data, Bioinformatics, Volume 40, Issue 3, March 2024, btae102, https://doi.org/10.1093/bioinformatics/btae102

Minimap2
Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574.

Samtools

Seurat

### **Contributors**
This package is developed and maintained by [Andre Sim](https://github.com/andredsim), [Min Hao Ling](https://github.com/lingminhao) and [Jonathan Goeke](https://github.com/jonathangoeke) at the Genome Institute of Singapore. If you want to contribute, please leave an issue. Thank you.
