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
We can test if this pipeline is installed correctly by running the following code using a small test set that comes with the Docker container.

``` 
nextflow run GoekeLab/bambu-SingleCell-Spatial \
  --reads reads_chr9_1_1000000.fastq.gz \
  --genome Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr9_1_1000000.fa \
  --annotation Homo_sapiens.GRCh38.91_chr9_1_1000000.gtf \
  --whitelist 737K-august-2016.txt.gz \ 
  --chemistry 10x5v2 --technology ONT --ncore 4 --outdir output \
  -with-singularity 
``` 

You can run the code above with your own dataset by replacing the arguments accordingly. The pipeline works on both the **long read sequencing** platform Oxford Nanopore Technologies (ONT) and PacBio at the single-cell and spatial level. The description for each argument is shown below: 

**1.reads** - a path to the read file (.fastq.gz) from a single-cell or spatial experiment run. The read name of this file has the follows the usual naming convention for reads from single-cell or spatial data: e.g. `ATCCGTCCAACGGGTA_TTGCTGGCGTGT#46f0ce76-6a12-4a12-a707-2ffef7be7594_+1of1`. The first 16 characters `ATCCGTCCAACGGGTA` refer to the cell barcode (CB) and the next 12 characters `TTGCTGGCGTGT` after the underscore refers to the unique molecular identifier (UMI).   

**2.genome** - a path to the genome fasta file. This should be the same file used for read alignment.

**3.annotation** - a path to the reference annotation file (.gtf), a txdb object, or annotations object prepared by prepareAnnotations(). When not provided, de novo transcript discovery is performed. (see [Bambu](https://github.com/GoekeLab/bambu) for more information)

**4.whitelist** - a path to the `gzipped` cell/spot barcode whitelist or the name of the barcode list ([single-cell](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-whitelist-) / [spatial](https://kb.10xgenomics.com/hc/en-us/articles/360041426992-Where-can-I-find-the-Space-Ranger-barcode-whitelist-and-their-coordinates-on-the-slide-)) provided by 10X Genomics. It is recommended to provide the correct whitelist as it is used as a reference to filter against low confidence putative barcodes, however you may also run this pipeline without providing it

<!--- | Choice of argument            | Name of barcodes list in 10x Genomics     | 
|-------------|:----------------| 
| visium-V1      |`visium-v1.txt` (The `V1` in `visium-V1 ` refers to the slide serial number for the Visium platform. The serial number may be replaced with the numbers that can range from 1 to 5) |
| 10x3v3 (default)      | `3M-february-2018.txt.gz`|
| 10x5v2      | `737-august-2016.txt`| ---> 

**5.technology** - the long read sequencing platform used to generate the reads [`ONT` or `PacBio`. default:`ONT`]

**6.chemistry** - type of 10X Genomics chemistry used in the single-cell or spatial experiment

| Choice of argument            | Description     | 
|:---|:----------------| 
| 10x3v2      | 10x 3' version 2|
| 10x3v3  (default)    | 10x 3' version 3|
| 10x5v2      | 10x 5' version 2|

**7.ncore** - the number of cores to use in the run (default: `12`)

**8.outdir** - the path to the output directory

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

### **Release History** 
TBD

### **Citation**
- bambuSC Paper
- Flexiplex Paper

### **Contributors**
This pipeline is developed and maintained by ... at the Genome Institute of Singapore. If you want to contribute, please leave an issue. Thank you.
