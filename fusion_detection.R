library(Biostrings)
library(data.table)
devtools::load_all("/mnt/software/bambu")

args <- commandArgs(trailingOnly = TRUE)
genome = args[[1]]
annotations = args[[2]]
jaffa_results = args[[3]]
sourceDir = args[[4]]
source(file.path(sourceDir,"fusion_detection_functions.R"))

cat('Load transcript sequence information')
geneSeq <- readDNAStringSet(file=genome)
listNames <- unlist(lapply(strsplit(names(geneSeq)," "),'[[',1))

jaffa_results.csv <- jaffa_results
jaffa_results_table <- fread(jaffa_results.csv)
jaffa_results_table$fusion.genes = jaffa_results_table$`fusion genes`
jaffa_results_table = jaffa_results_table[jaffa_results_table$classification == "HighConfidence",]
fusionGeneNames <- jaffa_results_table$fusion.genes


if(sum(gsub("chr", "", jaffa_results_table$chrom1) %in% listNames) > 
        sum(jaffa_results_table$chrom1 %in% listNames)){
    print("removing chr prefix from scaffolds")
    jaffa_results_table$chrom1 = gsub("chr", "", jaffa_results_table$chrom1)
    jaffa_results_table$chrom2 = gsub("chr", "", jaffa_results_table$chrom2)
}

anno.file <- annotations
anno.file <-prepareAnnotationsFromGTF(anno.file)

exonsByTx = anno.file
exonsByGene = split(unlist(exonsByTx), rep(mcols(exonsByTx)$GENEID, times = lengths(exonsByTx)))
exonsByGene <- endoapply(exonsByGene, unique)

ensemblAnnotations.genes = data.table(as_tibble(mcols(exonsByTx)) %>% 
                                         group_by(GENEID) %>% 
                                         summarise(ensembl_gene_id = GENEID[1],
                                                   hgnc_symbol = hgnc_symbol[1]))
ensemblAnnotations.txs = data.table(as_tibble(mcols(exonsByTx)) %>% rename(ensembl_gene_id = GENEID,
                                                                             ensembl_transcript_id = TXNAME))

fusionTx <- getFusionTxRange(fusionGeneNames, ensemblAnnotations.genes, exonsByGene, ensemblAnnotations.txs, exonsByTx,jaffa_results_table)
isNULL = sapply(fusionTx, is.null)
fusionAnnotations <- getFusionAnnotations(fusionGeneNames[!isNULL], anno.file,ensemblAnnotations.genes,ensemblAnnotations.txs,fusionTx[!isNULL])
full_annotations <- anno.file
filtered_annotations <- full_annotations[ensemblAnnotations.txs[hgnc_symbol %in% unlist(strsplit(fusionGeneNames, ":"))]$ensembl_transcript_id]
gr <- unlist(filtered_annotations)
df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr))

write.table(df, file="fusion.bed", quote=F, sep="\t", row.names=F, col.names=F)
writeToGTF(fusionAnnotations, file = "fusion.gtf")