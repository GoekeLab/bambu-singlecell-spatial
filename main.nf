#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.chemistry = "custom" //"10x3v3" "10x3v2" "10x5v2" "10x5v3"
params.technology = 'ONT' //"ONT" "PacBio"
params.whitelist = "NULL"
params.bambuPath = "bambu"
params.lowMemory = "FALSE"
params.ncore = 4
params.spatial = null

params.NDR = "NULL"
params.cleanReads = "TRUE"
params.keepChimericReads = "FALSE"
params.deduplicateUMIs = "TRUE"
params.bambuParams = tuple(params.cleanReads, params.keepChimericReads, params.deduplicateUMIs)
params.barcodeMap = "TRUE"
params.clusters = "auto"
params.resolution = 0.8
params.reads = null
params.bams = null
params.noEM = null

params.flexiplex_x = 'CTACACGACGCTCTTCCGATCT' //sequence Append flanking sequence to search for
params.flexiplex_b = '????????????????' //sequence Append the barcode pattern to search for
params.flexiplex_u = '????????????' //sequence Append the UMI pattern to search for
params.flexiplex_x2 = 'TTTTTTTTT' //sequence Append flanking sequence to search for
params.flexiplex_e = 1

params.outdir = "output" 

process flexiplex{ 

	publishDir "$params.outdir", mode: 'copy' 
	cpus params.ncore
	maxForks params.ncore
	
	input: 
	tuple val(sample), path(fastq), val(chemistry), val(technology), val(whitelist)

	output:
    tuple val(sample), path("new_reads.fastq.gz"), val(chemistry), val(technology), emit: fastq

	script:
	"""	
    chem="$chemistry"
    if [[ "$chemistry" == 10x3v3 ]]; then 
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT"
	elif [[ "$chemistry" == 10x3v2 ]]; then 
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTTTTTTT"
	elif [[ "$chemistry" == 10x5v2 ]]; then
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG"
    elif [[ "$chemistry" == 10x5v3 ]]; then
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTCTTATATGGG"
	fi

	gunzip -c $fastq > reads.fastq 
	flexiplex -p $params.ncore \$chem -f 0 reads.fastq
    if [[ $whitelist == "NULL" ]]; then
        python /mnt/software/main.py --outfile my_filtered_barcode_list.txt flexiplex_barcodes_counts.txt 
    else
        gunzip -c $whitelist > whitelist.txt 
	    python /mnt/software/main.py --outfile my_barcode_list.txt flexiplex_barcodes_counts.txt 
        awk '{print \$1}' whitelist.txt my_barcode_list.txt | sort | uniq -d > my_filtered_barcode_list.txt
	fi
    flexiplex -p $params.ncore -k my_filtered_barcode_list.txt \$chem -f 8 -e $params.flexiplex_e reads.fastq > new_reads.fastq
	gzip new_reads.fastq
    rm reads.fastq
    """
}

process minimap{ 
    publishDir "$params.outdir", mode: 'copy' 

	cpus params.ncore
	maxForks params.ncore
	
	input: 
	tuple val(sample),path(newfastq), val(chemistry), val(technology)
	path(genome)

	output: 
	tuple val(sample), path ('*.demultiplexed.bam') 

	script:
	""" 
	if [[ $technology == PacBio ]]; then 
		minimap2 -ax splice:hq -t $params.ncore -d ref.mmi $genome
		minimap2 -ax splice:hq -t $params.ncore -a ref.mmi $newfastq > demultiplexed.sam  
	
	else
		minimap2 -ax splice -k14 -t $params.ncore -d ref.mmi $genome
		minimap2 -ax splice -k14 -t $params.ncore -a ref.mmi $newfastq > demultiplexed.sam  
	fi

	samtools sort -@ $params.ncore demultiplexed.sam -o ${sample}.demultiplexed.bam 
	samtools index -@ $params.ncore ${sample}.demultiplexed.bam 

	rm demultiplexed.sam
	rm ref.mmi 
	"""
}

process bambu{ 
    publishDir "$params.outdir", mode: 'copy' 

	cpus params.ncore 
	maxForks 1
	
	input: 
	val(id)
    val(bam)
	path(genome)
	path(annotation)
    val(bambuPath)
    tuple val(cleanReads), val(keepChimericReads), val(deduplicateUMIs)
    val(NDR)
    val(barcode_map)
    val(whitelist)
    val(clusters)
    val(resolution)
    

	output: 
    tuple val ('combined'), path ('*readClassFile.rds'), path ('*quantData.rds')
	path ('*extended_annotations.rds') 
    path ('*.gtf')
    path ('*.mtx')
    path ('*.tsv')
    path ('*_clusters.rds'), emit: clusters

	script:
	""" 
	#!/usr/bin/env Rscript
    #.libPaths("/usr/local/lib/R/site-library")

    samples = "$bam"
    samples = gsub("[][]","", gsub(' ','', samples))
    samples = unlist(strsplit(samples, ','))

    ids = "$id"
    ids = gsub("[][]","", gsub(' ','', ids))
    ids = unlist(strsplit(ids, ','))
    if(length(ids)>1){runName = "combined_"
    }else{runName = paste0(ids, "_")}

    if(!isTRUE($barcode_map)){
        x = gsub("[][]","",gsub(' ','', "$barcode_map"))
        barcode_maps = unlist(strsplit(x, ','))
    } else {
        barcode_maps = TRUE
    }

    library(devtools)
    if("$bambuPath" == "bambu") {
        load_all("/mnt/software/bambu")
    } else {
        load_all("$bambuPath")
    }

    spatial = "$whitelist"
    if(spatial == "FALSE"){spatial = NULL}

	annotations <- prepareAnnotations("$annotation")

	# Transcript discovery and generate readGrgList for each cell
    readClassFile = bambu(reads = samples, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = barcode_maps, verbose = TRUE, assignDist = FALSE, lowMemory = as.logical("$params.lowMemory"), yieldSize = 10000000, sampleNames = ids, cleanReads = as.logical($cleanReads), dedupUMI = as.logical($deduplicateUMIs))
    saveRDS(readClassFile, paste0(runName, "_readClassFile.rds"))
    if(isFALSE($NDR)){
        extendedAnno = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    } else{
        extendedAnno = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, NDR = $NDR)
    }
    saveRDS(extendedAnno, paste0(runName, "_extended_annotations.rds"))
    rm(annotations)
    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE, spatial = spatial)
    saveRDS(se, paste0(runName, "_quantData.rds"))
    for(se.x in se){
        if(as.logical("$params.lowMemory")){
            writeBambuOutput(se.x, '.', prefix = metadata(se.x)\$sampleNames)
        } else{
            writeBambuOutput(se.x, '.', prefix = "combined_")
        }
        
    }
    #writeBambuOutput(do.call(cbind, se), '.')
    #write(runName, "runName.txt")

    #if no clustering provided, automatically cluster
    if("$clusters" == "auto"){
        clusters = list()
        cellMixs = list()
        source("${projectDir}/utilityFunctions.R")
        for(quantData in se){
            quantData.gene = transcriptToGeneExpression(quantData)
            for(sample in unique(colData(quantData)\$sampleName)){
                i = which(colData(quantData)\$sampleName == sample)
                counts = assays(quantData.gene)\$counts[,i]
                cellMix = clusterCells(counts, resolution = $resolution)
                x = setNames(names(cellMix@active.ident), cellMix@active.ident)
                names(x) = paste0(sample,"_",names(x))
                clusters = c(clusters, splitAsList(unname(x), names(x)))
                cellMixs = c(cellMixs, cellMix)
            }
        }
        saveRDS(cellMixs, paste0(runName, "_cellMixs.rds"))
    }
    if("$clusters" == "none"){
        clusters = NULL
    }
    saveRDS(clusters, paste0(runName, "_clusters.rds"))
	"""
}

process bambu_EM{

	publishDir "$params.outdir", mode: 'copy'

	cpus params.ncore
    maxForks 1

    input:
    tuple val(sample),path(readClassFile), path(quantData)
    path(extendedAnno)
    path(gtf)
    path(counts)
    path(metadata)
    path(clusters)
    path(genome)
    val(bambuPath)
    val(clusters)

    output:
    path ('*.rds')
    path ('*.mtx')
    path ('*.tsv')

    script:
    """
    #!/usr/bin/env Rscript
    #.libPaths("/usr/local/lib/R/site-library")
    library(devtools)
    if("$bambuPath" == "bambu") {
        load_all("/mnt/software/bambu")
    } else {
        load_all("$bambuPath")
    }
    if(".txt" %in% "$sample"){runName = readLines("$sample")
    } else{runName = "$sample"}
    
    extendedAnno <- readRDS("$extendedAnno")
    quantDatas = readRDS("$quantData")
    clusters = readRDS("$clusters")
    degBias = TRUE
    if(is.null(clusters)){degBias = FALSE}

    se = bambu( reads = "test.rds", 
                annotations = extendedAnno, 
                genome = "$genome", 
                quantData = quantDatas, 
                assignDist = FALSE, 
                ncore = $params.ncore, 
                discovery = FALSE, 
                quant = TRUE, 
                demultiplexed = TRUE, 
                verbose = FALSE, 
                opt.em = list(degradationBias = degBias), 
                clusters = clusters)

    saveRDS(se, paste0(runName, "_se.rds"))
    writeBambuOutput(se, path = ".", prefix = paste0(runName, "_EM_"),outputExtendedAnno = FALSE, 
        outputAll = FALSE, outputBambuModels = FALSE, outputNovelOnly = FALSE)
    """

}

process fusion_mode_JAFFAL{
    publishDir "$params.outdir", mode: 'copy'
    container ''

    input:
    tuple val(sample),path(fastq), val(chemistry), val(technology)
    val(jaffal_ref_dir)

    output:
    path('jaffa_results.csv')

    script:
    """
    fastaPath=\$(dirname $fastq)
    echo \$fastaPath
    fastaPath=\$(realpath \$fastaPath)
    fastq2=\$(realpath $fastq)
    apptainer run                             \
    -B $jaffal_ref_dir:/ref                      \
    -B \$fastaPath                      \
    docker://davidsongroup/jaffa:latest   \
    -n $params.ncore \
    /JAFFA/JAFFAL.groovy            \
    \$fastq2

    """
}

process fusion_mode_extract{
    publishDir "$params.outdir", mode: 'copy'

    input:
    path(jaffa_results)
    path(bam)
    path(genome)
    path(annotation)
    val(jaffal_ref_dir)

    output:
    path('fusionGene.fasta')
    path('sample1_fusion.bam')
    path('fusion.gtf')

    script:
    """
    Rscript ${projectDir}/fusion_detection.R $genome $annotation $jaffa_results ${projectDir}

    samtools view -bhL fusion.bed $bam | samtools fastq - > sample1_fusion.fastq
    minimap2 -ax splice -G2200k -N 5 --sam-hit-only -t $params.ncore fusionGene.fasta sample1_fusion.fastq | samtools sort -O bam -o sample1_fusion.bam -
    samtools index sample1_fusion.bam 
    """
}

process fusion_mode_bambu{
    publishDir "$params.outdir", mode: 'copy'

    input:
    path(fusionGeneScaffolds)
    path(bam)
    path(fusion_gtf)
    val(bambuPath)
    tuple val(cleanReads), val(keepChimericReads), val(deduplicateUMIs)

    output:
    path('*_fusion_readClassFile.rds')
    path('*_fusion_extended_annotations_NDR1.rds')
    path('*fusion_quantData.rds')

    script:
    """
    #!/usr/bin/env Rscript

    library(devtools)
    if("$bambuPath" == "bambu") {
        load_all("/mnt/software/bambu")
    } else {
        load_all("$bambuPath")
    }

	annotations <- prepareAnnotations("$fusion_gtf")

    readClassFile = bambu(reads = "$bam", annotations = annotations, genome = "$fusionGeneScaffolds", fusionMode = TRUE, ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = TRUE, assignDist = FALSE, lowMemory = as.logical("$params.lowMemory"), yieldSize = 10000000, cleanReads = as.logical($cleanReads), dedupUMI = as.logical($deduplicateUMIs))
    saveRDS(readClassFile, paste0("_fusion_readClassFile.rds"))
    extendedAnno_NDR1 = bambu(reads = readClassFile, annotations = annotations, genome = "$fusionGeneScaffolds", NDR = 1, fusionMode = TRUE, ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    saveRDS(extendedAnno_NDR1,paste0("_fusion_extended_annotations_NDR1.rds"))
    rm(annotations)
    se = bambu(reads = readClassFile, annotations = extendedAnno_NDR1, genome = "$fusionGeneScaffolds", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE)
    saveRDS(se, paste0("_fusion_quantData.rds"))
    for(se.x in se){
        writeBambuOutput(se.x, '.', prefix = paste0(metadata(se.x)\$sampleNames, "fusion"))
    }
    """
}

process fusion_mode_bambu_EM{

	publishDir "$params.outdir", mode: 'copy'

	cpus params.ncore
    maxForks 1

    input:
    path(readClassFile)
    path(extendedAnno)
    path(quantData)
    path(genome)
    val(bambuPath)
    path(clusters)

    output:
    path ('*.rds')
    path ('*.mtx')
    path ('*.tsv')

    script:
    """
    #!/usr/bin/env Rscript
    #.libPaths("/usr/local/lib/R/site-library")
    library(devtools)
    if("$bambuPath" == "bambu") {
        load_all("/mnt/software/bambu")
    } else {
        load_all("$bambuPath")
    }
    
    extendedAnno <- readRDS("$extendedAnno")
    quantDatas = readRDS("$quantData")
    clusters = readRDS("$clusters")
    degBias = TRUE
    if(is.null(clusters)){
        degBias = FALSE
    } else{
        sampleID = gsub("[AGCT]*\$", "", clusters[[1]][[1]][1], perl = TRUE)
        fusionID = paste0(colData(quantDatas[[1]])\$sampleName[1], "_")
        clusters[[1]] = gsub(sampleID, fusionID, clusters[[1]])
    }
    print(clusters)
    se = bambu( reads = "test.rds", 
                annotations = extendedAnno, 
                genome = "$genome", 
                quantData = quantDatas, 
                assignDist = FALSE, 
                ncore = $params.ncore, 
                discovery = FALSE, 
                quant = TRUE, 
                demultiplexed = TRUE, 
                verbose = FALSE, 
                opt.em = list(degradationBias = degBias), 
                clusters = clusters)

    saveRDS(se, paste0("_fusion_se.rds"))
    writeBambuOutput(se, path = ".", prefix = paste0("_fusion_EM_"), outputExtendedAnno = FALSE, 
        outputAll = FALSE, outputBambuModels = FALSE, outputNovelOnly = FALSE)
    """

}

// This is the workflow to execute the process 
workflow {    
	if (params.reads) {
        //User can provide either 1 .fastq file or a .csv with .fastq files
        fastq = file(params.reads, checkIfExists:true)
        if(fastq.getExtension() == "csv") {
            //TODO NEED TO CHECK THAT ALL CHEMISTRY AND TECHNOLOGY VALUES ARE VALID
            //TODO NEED TO INCLUDE WHITELIST

            readsChannel = Channel.fromPath(params.reads) \
                    | splitCsv(header:true, sep:',') \
                    | map { row-> tuple(row.sample, file(row.fastq), 
                                        row.containsKey("chemistry") ? row.chemistry : params.chemistry,
                                        row.containsKey("technology") ? row.technology : params.technology,
                                        row.containsKey("whitelist") ? row.whitelist : params.whitelist,
                                        row.containsKey("barcode_map") ? row.barcode_map :  params.barcodeMap,
                                        row.containsKey("clusters") ? row.whitelist : params.clusters) }        
            flexiplex_out_ch = flexiplex(readsChannel.map{it[0..4]})
            minimap_out_ch = minimap(flexiplex_out_ch, "$params.genome")
            barcodeMaps = readsChannel.collect{it[6]}
            barcodeMaps2 = barcodeMaps.map { it == null ? it : params.barcodeMap }
            whitelists = readsChannel.collect{it[5]}
            whiteLists2 = whitelists.map { it == null ? it : params.whitelist }
            clusters = readsChannel.collect{it[7]}
            clusters2 = clusters.map { it == null ? it : params.clusters }
        }
        else {  
            //if (params.chemistry != "10x3v3"| params.chemistry != "10x3v2" | params.chemistry != "10x5v2"){ exit 1, "--chemistry must be one of (3prime-v3/3prime-v2/5prime-v2)" }
            //if (params.technology == false) { exit 1, "--technology must be one of (ONT/PacBio)" }
            readsChannel = Channel.fromPath(params.reads)
            readsChannel = readsChannel
                .map {tuple("Bambu", it, params.chemistry, params.technology, params.whitelist)}
            flexiplex_out_ch = flexiplex(readsChannel)
            minimap_out_ch = minimap(flexiplex_out_ch, "$params.genome")
            barcodeMaps2 = params.barcodeMap
            whiteLists2 = params.whitelist
            clusters2 = params.clusters

        }
        sampleIds = minimap_out_ch.collect{it[0]}
        bamsFiles = minimap_out_ch.collect{it[1]}
        if (params.fusionMode) {
            fusion_mode_JAFFAL_out_ch = fusion_mode_JAFFAL(readsChannel.map{it[0..3]}, "$params.jaffal_ref_dir")
            fusion_mode_extract_out_ch = fusion_mode_extract(fusion_mode_JAFFAL_out_ch, bamsFiles.flatten(), "$params.genome", "$params.annotation", "$params.jaffal_ref_dir")
            fusion_mode_bambu_out_ch = fusion_mode_bambu(fusion_mode_extract_out_ch, "$params.bambuPath", params.bambuParams)
        }
     }
    else{ //When starting from bam
        bam = file(params.bams, checkIfExists:true)
        if(bam.getExtension() == "csv") {
            bamsChannel = Channel.fromPath(params.bams) \
                    | splitCsv(header:true, sep:',') \
                    | map { row-> 
                                def barcodeMap = row.containsKey("barcode_map") ? row.barcode_map : null
                                def spatialWhitelist = row.containsKey("spatial_whitelist") ? row.spatial_whitelist : null
                                def clusters = row.containsKey("clusters") ? row.clusters : null
                                tuple(row.sample, file(row.bam), barcodeMap, spatialWhitelist, clusters) }
            barcodeMaps = bamsChannel.collect{it[2]}
            whiteLists = bamsChannel.collect{it[3]}
            clusters = bamsChannel.collect{it[4]}

            barcodeMaps2 = barcodeMaps.map { it == null ? it : params.barcodeMap }
            clusters2 = clusters.map { it == null ? it : params.clusters }
            whiteLists2 = whiteLists.map { it == null ? it : params.whitelist }
        }
        if(bam.getExtension() == "bam"){
            bamsChannel = Channel.fromPath(params.bams)
            bamsChannel = bamsChannel
                    .map {["Bambu", it]}

            barcodeMaps2 = params.barcodeMap
            whiteLists2 = params.whitelist
            clusters2 = params.clusters
        }
        sampleIds = bamsChannel.collect{it[0]}
        bamsFiles = bamsChannel.collect{it[1]}



    }
    if(!params.spatial){
        whiteLists2 = "FALSE" 
    }
    bambu_out_ch = bambu(sampleIds, bamsFiles, "$params.genome", "$params.annotation", "$params.bambuPath", params.bambuParams,"$params.NDR",barcodeMaps2, whiteLists2, clusters2, "$params.resolution")
    if(!params.noEM){
        bambuEM_out_ch = bambu_EM(bambu_out_ch, "$params.genome", "$params.bambuPath", bambu_out_ch.clusters)
    }
    if (params.fusionMode && !params.noEM && params.reads) {
        fusion_mode_bambu_EM(fusion_mode_bambu_out_ch, "$params.genome", "$params.bambuPath", bambu_out_ch.clusters)
    }
}

