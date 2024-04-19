# **Context-Aware Transcript Quantification from Long Read Single-Cell and Spatial Transcriptomics data**
This is a pipeline developed for context-aware transcript quantification from long read single-cell and spatial transcriptomics data. The pipeline consists of two steps: first we use [flexiplex](https://davidsongroup.github.io/flexiplex/) for cell or spatial barcode discovery and demultiplexing. Then we use  [Bambu](https://github.com/GoekeLab/bambu/tree/BambuDev) that has been adapted for single-cell or spatial data for transcript discovery and quantification. For user's convenience, we wrapped these tools using Nextflow running on a Docker container so that they can be executed using an one-line code. The output from this pipeline can then be used for visualization and downstream analysis such as differential isoform expression or isoform switching event detection.

### **Content** 
- [Installation](#installation)
- [General Usage](#General-Usage)
- [Release History](#Release-History)
- [Citation](#Citation)
- [Contributors](#Contributors)


### **Installation** 
To run this pipeline, you will need Nextflow and Docker (or Singularity if you are working on a high performance computing clusters). Latest version is recommended. 

### **General Usage** 
You can test if this pipeline is installed correctly by running the following code using a small test set that comes with the Docker container.

``` 
nextflow run GoekeLab/bambu-SingleCell-Spatial \
  --reads reads_chr9_1_1000000.fastq.gz \
  --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \
  --whitelist 737K-august-2016.txt.gz \ 
  --chemistry 10x5v2 --ncore 4 --outdir output \
  -with-singularity 
``` 

You can run the code above with your own dataset by replacing the arguments accordingly. The pipeline works on both the **long read sequencing** platform Oxford Nanopore Technologies (ONT) and PacBio at the single-cell and spatial level. The description for each argument is shown below: 

### **Required arguments**

The BambuSC pipeline can be started from from raw reads (fastq.gz) or from a demultiplexed .bam file if you have already produced these (from earlier runs of this pipeline or other upstream tools). Therefore either the --reads or --bams argument is mandatory depending on your input files.

**--reads** [PATH] - a path to one read file (.fastq.gz) from a single-cell or spatial experiment run. 

Multiple Samples:

If running multiple samples you can instead path to a csv sample sheet containing the parameters neede for each sample. Samples with the same 'sample' name will be combined together and treated as one sample (for example when combining technical replicates that have the same barcodes). See the following argument descriptions on how the samplesheet can be filled out. If the chemistry, technology, or whitelist columns are not provided, or if there is a missing entry, the pipeline will use the input from the --chemistry --technology and --whitelist arguments for the samples. 

Multiple Sample Example:
| sample            | fastq     | chemistry      | technology | whitelist |
|:---|:----------|:------|:------|:------| 
| condition1      | path/to/1.fastq.gz | 10x3v2 | ONT | path/to/whitelist.gz |
| uniqueSample1      | path/to/2.fastq.gz | -x CTACACGACGCTCTTCCGATCT -b ???????????? -u ???????????? -x TTTTTTTTT | PacBio | path/to/whitelist2.gz |

**--bams** [PATH] - If you did demultiplexing and alignment seperately, you can start the pipeline at the bambu step. Requires a path to the alignment file (.bam). The bam file must be demultiplexed in one of the following 3 ways:

1. The read names have the barcode and UMI sequences appended to the start of the name: e.g. `ATCCGTCCAACGGGTA_TTGCTGGCGTGT#46f0ce76-6a12-4a12-a707-2ffef7be7594_+1of1`. The first set of characters (`ATCCGTCCAACGGGTA`) followed by an '_' refer to the cell barcode (CB) and the next set of characters after the underscore  (`TTGCTGGCGTGT`) followed by a '#' refers to the unique molecular identifier (UMI). If there is no UMI in your data this set can be empty and no UMI deduplication will be performed, but ensure the underscore and hash are still present so that the bam file name can be correctly parsed.

2. The barcode is located in the BC tag in the bam file, and the UMI is located in the UG tag

3. The --barcode_map argument is provided in (or as a column in the samplesheet) to a .tsv or .csv with three columns. All reads not present in the table will be discarded from the analysis. The UMI column is optional if you do not need UMI deduplication or do not have UMIs in your data.
Example:

| Read Name            | Barcode     | UMI      |
|:---|:----------|:------| 
| #46f0ce76-6a12-4a12-a707-2ffef7be7594      | ATCCGTCCAACGGGTA | TTGCTGGCGTGT |

Multiple Samples:

If running multiple samples you can instead path to a csv sample sheet containing the parameters neede for each sample. Samples with the same 'sample' name will be combined together and treated as one sample (for example when combining technical replicates that have the same barcodes). See the following argument descriptions on how the samplesheet can be filled out. Providing a barcode map is optional, and only required if the bam file does not contain the barcode information. Providing a spatial_whitelist is only required if using spatial data, and each sample uses a different whitelist/coordinates, otherwise it will use the whitelist providied by the --whitelist argument.

Multiple Sample Example:
| sample            | bam     | barcode_map | spatial_whitelist
|:---|:----------|:------|:------|  
| condition1      | path/to/demultiplexed.bam | | |
| condition2      | path/to/not_demultiplexed.bam | path/to/condition2_barcode_map.csv | |

**--genome** [PATH] - a path to the genome fasta file. in the event you are running the pipeline using a bam file, this should be the same file used for read alignment.

**--annotation** [PATH] - a path to the reference annotation file (.gtf), a txdb object, or annotations object prepared by prepareAnnotations(). When not provided, de novo transcript discovery is performed. (see [Bambu](https://github.com/GoekeLab/bambu) for more information)

<!--- | Choice of argument            | Name of barcodes list in 10x Genomics     | 
|-------------|:----------------| 
| visium-V1      |`visium-v1.txt` (The `V1` in `visium-V1 ` refers to the slide serial number for the Visium platform. The serial number may be replaced with the numbers that can range from 1 to 5) |
| 10x3v3 (default)      | `3M-february-2018.txt.gz`|
| 10x5v2      | `737-august-2016.txt`| ---> 

**--chemistry** [10x3v2/10x3v3/10x5v2/STRING] - type of 10X Genomics chemistry used in the single-cell or spatial experiment

| Choice of argument            | Description     | 
|:---|:----------------| 
| 10x3v2      | 10x 3' version 2|
| 10x3v3      | 10x 3' version 3|
| 10x5v2      | 10x 5' version 2|

Advanced use: If your sample was not barcoded by one of the above chemistries you can provide a custom search string as detailed in the flexiplex documentation (https://davidsongroup.github.io/flexiplex/). e.g --chemistry '-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTTTTTTT'


### **Optional arguments**
**--flexiplex_e** [INT] - the edit distance used by flexiplex for barcode identification (default: 1)

**--technology** [ONT/PacBio]- the long read sequencing platform used to generate the reads [`ONT` or `PacBio`. default:`ONT`]
This changes the stringency of the alignment based on minimap2 recommendations.

**--NDR** [FLOAT] - the NDR threshold used for Bambu transcript discovery. A value between 0 - 1, with lower values being more precise but less sensitive. Typically a value of 0.1 is suitable for most analyses. (default: automatically determined based on data)

**--BambuDiscoveryParameters** [STRING] - A comma seperated sting containing optional parameters for bambu's transcript discovery opt.discovery argument e.g --BambuDiscoveryParameters 'remove.subsetTx = FALSE, min.readFractionByGene = 0'

**--BambuQuantifyParameters** [STRING] - A comma seperated sting containing optional parameters for bambu's transcript quantification opt.em argument e.g --BambuQuantifyParameters 'maxiter = 20000, conv = 0.0002'

**--whitelist** [PATH] - a path to the `gzipped` cell/spot barcode whitelist or the name of the barcode list ([single-cell](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-) / [spatial](https://kb.10xgenomics.com/hc/en-us/articles/360041426992-Where-can-I-find-the-Space-Ranger-barcode-whitelist-and-their-coordinates-on-the-slide-)) provided by 10X Genomics. It is recommended to provide the correct whitelist as it is used as a reference to filter against low confidence putative barcodes, however you may also run this pipeline without providing it

**--ncore** [INT] -  the number of cores to use in the run (default: `12`)

**--keepChimericReads** [TRUE/FALSE] - determines if chimeric reads detected by flexiplex will be kept. Setting this to false will reduce the number of reads used in the analysis, but speed up the pipeline (default: 'true')

**--cleanReads [TRUE/FALSE] - determines if (default: 'true')

**--outdir** [PATH] - the path to the output directory

**--bambuPath** [PATH] - For advanced users. The path to a local version of Bambu which will be used

After the run, all the outputs will be stored in the folder specified by the `--outdir` parameter. The output files are described below: 

| Output file name                | Description                                                             |
|:----------------------|:------------------------------------------|
| extended_annotations.gtf        | Extended transcript & gene annotations for the genome using long reads data.        |
| sparse_counts_transcript.mtx.gz           | Total read counts estimates for each transcript in each sample (sparse matrix format).        |
| sparse_CPM_transcript.mtx.gz              | Counts per million (CPM) estimates for each transcript in each sample (sparse matrix format). |
| sparse_fullLengthCounts_transcript.mtx.gz | Full length read counts estimates for each transcript in each sample (sparse matrix format).  |
| sparse_uniqueCounts_transcript.mtx.gz                | Unique read counts estimates for each transcript in each sample (sparse matrix format).       |
| sparse_counts_gene.mtx.gz                 | Gene read counts estimates for each transcript in each sample.         |
| txANDGenes.tsv.gz                 | Gene ID associated to each transcript arranged as in the transcript count estimates          |
| genes.tsv.gz                 | Gene ID arranged as in the gene count estimates          |

You may then use these count matrices for downstream analysis using tools like [Seurat](https://satijalab.org/seurat/) or [scanpy](https://www.google.com/search?q=scanpy&oq=scanpy&aqs=chrome..69i57.866j0j7&sourceid=chrome&ie=UTF-8). The following short tutorials demonstrate how to analyse the long read [single-cell]() and [spatial]() data using [Seurat](https://satijalab.org/seurat/). 

### **Multiple Sample Analysis** ###

If you would like to analysis multiple samples or replicates, please note the following differences in the arguments. 
The input to the --reads/--bams arguments needs to be a .csv file with all columns present. See the arguments section to ensure the formatting of the samplesheet.csv is correct and contains all the required information.
The sample column is used to distinguish different replicates and the sample name will be appended to the barcodes from that sample. Therefore if you would like to combine samples that share barcodes you can assign them the same sample name and their counts will be comined. 
The same --genome  and --annotation argument inputs will be applied to all samples. 

### **Spatial Analysis** ##

To run spatial data using Bambu, the X and Y coordinate mapping for each barcode to a spot needs to be provided as two additional columns in the whitelist (tab seperated). This file can generally be sourced from the provider of the library preperation

Example (**Headers are descriptive only, there should be no headers in the input file!**):
| Barcode            | X coordinate | Y coordinate| 
|:---|:---|:---|
| AAACAACGAATAGTTC | 17 | 1 |
| AAACAAGTATCTCCCA  | 103 | 51 |
| AAACAATCTACTAGCA | 44 | 4 |

### **Fusion Transcript Analysis** ##
TBD

### **Additional Information**

UMI correction is done at the barcode level. The longest read for each unique barcode-UMI combination is kept for analysis.

To reduce the impact of chimeric reads, only the first alignment after the barcode will be used for analysis.

### **Release History** 
TBD

### **Citation**
- bambuSC Paper
- Flexiplex Paper

### **Contributors**
This pipeline is developed and maintained by ... at the Genome Institute of Singapore. If you want to contribute, notice bugs, or have suggestions, please leave a GitHub issue. Thank you.
