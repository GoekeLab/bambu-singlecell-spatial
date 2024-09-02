prepareAnnotationsFromGTF <- function(file) {
    if (missing(file)) {
        stop("A GTF file is required.")
    } else {
        data <- utils::read.delim(file, header = FALSE, comment.char = "#")
        colnames(data) <- c("seqname", "source", "type", "start", "end",
            "score", "strand", "frame", "attribute")
        data <- data[data$type == "exon", ]
        data$strand[data$strand == "."] <- "*"
        data$GENEID <- gsub("gene_id (.*?);.*", "\\1", data$attribute)
        data$TXNAME <- gsub(".*transcript_id (.*?);.*", "\\1", data$attribute)
        data$hgnc_symbol <- gsub(".*gene_name (.*?);.*","\\1",data$attribute)
        multiTxCheck <- as_tibble(data) %>% select(seqname, GENEID) %>% distinct() %>% group_by(GENEID) %>% 
            mutate(n=n(), id=paste0('-',row_number()))
        if(any(multiTxCheck$n>1)) { # identical TXNAMES
            warning('Annotations contain duplicated transcript names
                        Transcript names will be made unique')
            uniqueNamesTbl <- as_tibble(data) %>% 
                select(seqname, TXNAME, GENEID) %>% 
                mutate(order=row_number()) %>% left_join(multiTxCheck) %>%
               mutate(gene_unique = paste0(GENEID, ifelse(n==1,'', id)),
                      tx_unique = paste0(TXNAME, ifelse(n==1,'', id))) %>%
                arrange(order)
            
            data$TXNAME_Original <- data$TXNAME
            data$GENEID_Original <- data$GENEID
            data$TXNAME <- uniqueNamesTbl$tx_unique
            data$GENEID <- uniqueNamesTbl$gene_unique
            }
        geneData <- unique(data[, c("TXNAME", "GENEID","hgnc_symbol")])
        grlist <- makeGRangesListFromDataFrame(
        data[, c("seqname", "start", "end", "strand", "TXNAME")],
            split.field = "TXNAME", keep.extra.columns = TRUE)
        grlist <- grlist[IRanges::order(start(grlist))]
        unlistedExons <- unlist(grlist, use.names = FALSE)
        partitioning <- PartitioningByEnd(cumsum(elementNROWS(grlist)),
            names = NULL)
        txIdForReorder <- togroup(PartitioningByWidth(grlist))
        exon_rank <- lapply(elementNROWS(grlist), seq, from = 1)
        exon_rank[which(unlist(unique(strand(grlist))) == "-")] <- lapply(
            exon_rank[which(unlist(unique(strand(grlist))) == "-")], rev
            ) # * assumes positive for exon ranking
        names(exon_rank) <- NULL
        unlistedExons$exon_rank <- unlist(exon_rank)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            unlistedExons$exon_rank)]
        # exonsByTx is always sorted by exon rank, not by strand,
        # make sure that this is the case here
        unlistedExons$exon_endRank <- unlist(lapply(elementNROWS(grlist),
            seq, to = 1), use.names = FALSE)
        unlistedExons <- unlistedExons[order(txIdForReorder,
            start(unlistedExons))]
        mcols(unlistedExons) <- mcols(unlistedExons)[, c("exon_rank",
            "exon_endRank")]
        grlist <- relist(unlistedExons, partitioning)
        # sort the grlist by start position, ranked by exon number
        mcols(grlist) <- DataFrame(geneData[(match(names(grlist),
            geneData$TXNAME)), ])
        mcols(grlist)$txid <- seq_along(grlist)
        minEqClasses <- getMinimumEqClassByTx(grlist)
        if(!identical(names(grlist),minEqClasses$queryTxId)) warning('eq classes might be incorrect')
        mcols(grlist)$eqClassById <- minEqClasses$eqClassById
    }
    return(grlist)
}

getFusionAnnotations <- function(fusionGeneNames, anno.file,ensemblAnnotations.genes,ensemblAnnotations.txs,fusionTx){
    annotationRanges <- anno.file#prepareAnnotations(anno.file)
    fusionAnnotation <- lapply(fusionGeneNames, function(s){
    print(s)
    genevec <- unlist(strsplit(s, ":"))
    valid <- do.call("c",lapply(genevec, function(g){
        geneid <- ensemblAnnotations.genes[hgnc_symbol==g]$ensembl_gene_id
        if(length(geneid)>1){ return(FALSE)}
        else{return(TRUE)}
    }))
    if(!all(valid)){
      print("multiple gene ids for hgnc symbol, skipping")
      print(s)
      return(NULL)
    }
    tmp_range <- do.call("c",lapply(genevec, function(g){
    geneid <- ensemblAnnotations.genes[hgnc_symbol==g]$ensembl_gene_id
    txid <- ensemblAnnotations.txs[ensembl_gene_id==geneid]$ensembl_transcript_id
    print(length(txid))
    tmp <- do.call("c",lapply(txid, function(t){
      print(t)
      tmp <- annotationRanges[t]
      fusionTmp <- fusionTx[[s]][[t]]
      seqlevels(tmp, pruning.mode = "tidy") <- as.character(unique(seqnames(tmp)))
      seqlevels(tmp) <- as.character(unique(seqnames(fusionTmp)))                                 
      start(tmp[[1]]) <- start(fusionTmp)
      end(tmp[[1]]) <- end(fusionTmp)
      # if negative strand
      if(unique(strand(tmp[[1]]))=="-"){
        tmp[[1]]$exon_rank <- max(tmp[[1]]$exon_rank)-tmp[[1]]$exon_rank+1
        tmp[[1]]$exon_endRank <- max(tmp[[1]]$exon_endRank)-tmp[[1]]$exon_endRank+1
      }
      strand(tmp[[1]]) <- strand(fusionTmp)
      return(tmp)
    }))
    return(tmp)
  }))
  return(tmp_range)
})
keep.id <- which(!sapply(fusionAnnotation, is.null))
fusionAnnotation <- fusionAnnotation[keep.id]
fusionAnnotation <- do.call("c", fusionAnnotation)
return(fusionAnnotation)
}

getFusionTxRange <- function(fusionGeneNames, ensemblAnnotations.genes, exonsByGene, ensemblAnnotations.txs, exonsByTx, jaffa_results_table){
    fusionTx <- lapply(seq_along(fusionGeneNames), function(s){
  print(s)
  tmp <- fusionGeneNames[s]
  genevec <- unlist(strsplit(tmp, ":"))
  if(any(!(genevec %in% ensemblAnnotations.genes$hgnc_symbol))){ return(NULL) }
  geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[1]]$ensembl_gene_id
  tmp_range <- exonsByGene[geneid][[1]]
  min_start <- min(start(tmp_range)) 
  length_g1 <- max(end(tmp_range)) - min_start + 1
  tmp_range <- lapply(seq_along(genevec), function(g){
    geneid <- ensemblAnnotations.genes[hgnc_symbol==genevec[g]]$ensembl_gene_id
    tmp_range <- exonsByGene[geneid][[1]]
    min_start <- min(start(tmp_range)) 
    length_g_tmp <- max(end(tmp_range)) 
    strand.info <- as.character(unique(strand(tmp_range)))
    txid <- ensemblAnnotations.txs[(ensembl_gene_id == geneid)]$ensembl_transcript_id
    new_range <- lapply(txid, function(t){
      tmp_range <- exonsByTx[t][[1]]
      if(g==1){
        scoreNum = ifelse(strand.info=="-", length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base1+1,
                          jaffa_results_table[fusion.genes == tmp]$base1-min_start+1)
      }else{
        scoreNum = ifelse(strand.info=="-",length_g_tmp-jaffa_results_table[fusion.genes == tmp]$base2+1+length_g1,
                          jaffa_results_table[fusion.genes == tmp]$base2-min_start+1+length_g1)
      }
      if(strand.info=="-"){
        new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                             ranges = IRanges(start = length_g_tmp - end(tmp_range) + 1, 
                                              end = length_g_tmp - start(tmp_range) + 1 ),
                             strand = Rle("+",length(tmp_range)),
                             score = scoreNum)
      }else{
        new_range <- GRanges(seqnames=Rle(tmp,length(tmp_range)),
                             ranges = IRanges(start = start(tmp_range) - min_start+1, 
                                              end = end(tmp_range) - min_start+1),
                             strand = Rle(strand.info,length(tmp_range)),
                             score = scoreNum)
      }
      if(g == 2){
        new_range <- GenomicRanges::shift(new_range, shift = length_g1)
        
      }
      return(new_range)  
    })
    names(new_range) <- txid
    return(new_range)
  })
  
  return(do.call("c",tmp_range))
})
names(fusionTx) <- fusionGeneNames
return(fusionTx)
}
