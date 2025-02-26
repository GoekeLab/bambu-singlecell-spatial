library(Biostrings)
library(data.table)
devtools::load_all("/mnt/software/bambu")

args <- commandArgs(trailingOnly = TRUE)
genome = args[[1]]
annotations = args[[2]]
jaffa_results = args[[3]]
# sourceDir = args[[4]]
# source(file.path(sourceDir,"fusion_detection_functions.R"))
source(args[[4]])

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

fusionGeneNames.unique = unique(fusionGeneNames)

fusionGeneNames.toPrint = unlist(sapply(seq_along(fusionGeneNames.unique), function(s){
  tmp <- fusionGeneNames.unique[s]
  genevec <- unlist(strsplit(tmp, ":"))
  if(all(genevec %in% ensemblAnnotations.genes$hgnc_symbol)){
    return(tmp)
  }}))

fusionGeneSequence <- unlist(lapply(seq_along(fusionGeneNames.toPrint), function(s){
    tmp <- fusionGeneNames.toPrint[s]
    genevec <- unlist(strsplit(tmp, ":"))
    paste(unlist(lapply(seq_along(genevec), function(g){
        geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[g]]$ensembl_gene_id
        tmp_range <- exonsByGene[geneid][[1]]
        seq_pos <- min(start(tmp_range)):max(end(tmp_range))
        seqChar <- geneSeq[[match(as.character(unique(seqnames(tmp_range))), listNames)]][seq_pos]
        if(as.character(unique(strand(tmp_range))) == "-"){
            seqChar <- reverseComplement(seqChar)
        }
        as.character(seqChar)
    })), collapse = "")    
}))


sink("fusionGene.fasta")
noprint <- lapply(seq_along(fusionGeneNames.toPrint), function(s){
  cat(paste0(">",fusionGeneNames.toPrint[s]), " \n")
  cat(fusionGeneSequence[s]," \n")
})
sink()

# scaffoldNames = names(geneSeq)
# scaffoldNames = gsub("\\s.*", "", scaffoldNames, perl = TRUE)
# geneSeq2 = geneSeq[match(jaffa_results_table$contig, scaffoldNames)]
# sink("fusionGene.fasta")
# noprint <- lapply(seq_along(geneSeq2), function(s){
#   cat(paste0(">",jaffa_results_table$`fusion genes`[s]), " \n")
#   cat(as.character(geneSeq2[s])," \n")
# })
# sink()