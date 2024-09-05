#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.chemistry = "custom" //"10x3v3" "10x3v2" "10x5v2"
params.technology = 'ONT' //"ONT" "PacBio"
params.whitelist = "NULL"
params.bambuPath = "bambu"
params.lowMemory = "FALSE"
params.ncore = 4
params.spatial = "FALSE"

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
    tuple val(sample), path("new_reads.fastq"), val(chemistry), val(technology), emit: fastq

	script:
	"""	
    chem="$chemistry"
    if [[ $chemistry == 10x3v3 ]]; then 
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT"
	elif [[ $chemistry == 10x3v2 ]]; then 
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTTTTTTT"
	elif [[ $chemistry == 10x5v2 ]]; then
        chem="-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG"
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
    val(lowMemory)
    

	output: 
    tuple val ('combined'), path ('*readClassFile.rds'), path ('*quantData.rds')
	path ('*extended_annotations.rds') 
    path ('*extended_annotations_NDR1.rds') 
    path ('*extended_annotations_NDR1.gtf') 
    path ('*.mtx')
    path ('*.tsv')

	script:
	""" 
	#!/usr/bin/env Rscript

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
    
    
    print(samples)
    print(runName)
    print(barcode_maps)
    library(devtools)
    if("$bambuPath" == "bambu") {
        load_all("/mnt/software/bambu")
    } else {
        
        load_all("$bambuPath")
    }

	annotations <- prepareAnnotations("$annotation")

	# Transcript discovery and generate readGrgList for each cell
    readClassFile = bambu(reads = samples, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = barcode_maps, verbose = TRUE, assignDist = FALSE, lowMemory = $params.lowMemory, yieldSize = 10000000, sampleNames = ids, cleanReads = as.logical($cleanReads), dedupUMI = as.logical($deduplicateUMIs))
    saveRDS(readClassFile, paste0(runName, "_readClassFile.rds"))
    if(isFALSE($NDR)){
        extendedAnno = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    } else{
        extendedAnno = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, NDR = $NDR)
    }
    saveRDS(extendedAnno, paste0(runName, "_extended_annotations.rds"))
    extendedAnno_NDR1 = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", NDR = 1, ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    saveRDS(extendedAnno_NDR1,paste0(runName, "_extended_annotations_NDR1.rds"))
    writeToGTF(extendedAnno_NDR1, paste0(runName, "extended_annotations_NDR1.gtf"))
    rm(extendedAnno_NDR1)
    rm(annotations)
    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE, spatial = $whitelist)
    saveRDS(se, paste0(runName, "_quantData.rds"))
    writeBambuOutput(do.call(cbind, se), '.')
    #write(runName, "runName.txt")
	"""
}

process bambu_EM{

	publishDir "$params.outdir", mode: 'copy'

	cpus params.ncore
    maxForks 1

    input:
    tuple val(sample),path(readClassFile), path(quantData)
    path(extendedAnno)
    path(extendedAnno_NDR1)
    path(extended_annotations_NDR1_gtf) 
    path(counts)
    path(metadata)
    path(genome)
    val(bambuPath)
    val(clusters)
    val(resolution)

    output:
    path ('*.rds')
    path ('*.mtx')
    path ('*.gtf')
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

    #if no clustering provided, automatically cluster
    if("$clusters" == "auto"){
        clusters = list()
        cellMixs = list()
        source("${projectDir}/utilityFunctions.R")
        for(quantData in quantDatas){
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
        saveRDS(clusters, paste0(runName, "_clusters.rds"))
        saveRDS(cellMixs, paste0(runName, "_cellMixs.rds"))
    }
    if("$clusters" == "none"){
        clusters = NULL
    }

    se = bambu( reads = "test.rds", 
                annotations = extendedAnno, 
                genome = "$genome", 
                quantData = quantData, 
                assignDist = FALSE, 
                ncore = $params.ncore, 
                discovery = FALSE, 
                quant = TRUE, 
                demultiplexed = TRUE, 
                verbose = FALSE, 
                opt.em = list(degradationBias = FALSE), 
                clusters = clusters)

    saveRDS(se, paste0(runName, "_se.rds"))
    writeBambuOutput(se, path = ".", prefix = paste0(runName, "_EM_"))
    """

}

process fusion_mode{
    input:

    output:

    script:
    """
    bpipe run -p refBase=$jaffal_ref_dir $jaffal_ref_dir/JAFFAL.groovy $fastq
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
                                        row.containsKey("barcode_map") ? row.whitelist : params.whitelist) }
            flexiplex_out_ch = flexiplex(readsChannel)
            minimap_out_ch = minimap(flexiplex_out_ch, "$params.genome")
            barcodeMaps = readsChannel.collect{it[5]}
            barcodeMaps2 = barcodeMaps.map { it == null ? it : "TRUE" }
            whiteLists2 = params.whitelist
            //params.bams = minimap_out_ch
        }
        else {  
            //if (params.chemistry != "10x3v3"| params.chemistry != "10x3v2" | params.chemistry != "10x5v2"){ exit 1, "--chemistry must be one of (3prime-v3/3prime-v2/5prime-v2)" }
            //if (params.technology == false) { exit 1, "--technology must be one of (ONT/PacBio)" }
            readsChannel = Channel.fromPath(params.reads)
            readsChannel.set{ch_test}
            ch_test = ch_test
                .map {tuple("Bambu", it, params.chemistry, params.technology, params.whitelist)}
            flexiplex_out_ch = flexiplex(ch_test)
            minimap_out_ch = minimap(flexiplex_out_ch, "$params.genome")
            barcodeMaps2 = params.barcodeMap
            whiteLists2 = params.whitelist
            //params.bams = minimap_out_ch
        }
        sampleIds = minimap_out_ch.collect{it[0]}
        bamsFiles = minimap_out_ch.collect{it[1]}
     }
    else{ //When starting from bam
        bam = file(params.bams, checkIfExists:true)
        if(bam.getExtension() == "csv") {
            bamsChannel = Channel.fromPath(params.bams) \
                    | splitCsv(header:true, sep:',') \
                    | map { row-> 
                                def barcodeMap = row.containsKey("barcode_map") ? row.barcode_map : null
                                def spatialWhitelist = row.containsKey("spatial_whitelist") ? row.spatial_whitelist : null
                                tuple(row.sample, file(row.bam), barcodeMap, spatialWhitelist) }
        }
        if(bam.getExtension() == "bam"){
            bamsChannel = Channel.fromPath(params.bams)
            bamsChannel = bamsChannel
                    .map {["Bambu", it]}
        }
        sampleIds = bamsChannel.collect{it[0]}
        bamsFiles = bamsChannel.collect{it[1]}
        barcodeMaps = bamsChannel.collect{it[2]}
        whiteLists = bamsChannel.collect{it[3]}

        barcodeMaps2 = barcodeMaps.map { it == null ? it : "TRUE" }

        if(params.whitelist == "null" & whiteLists.any { it == '' }){
            echo "Missing whitelist entries in samplesheet or no --whitelist provided"
            exit 1
        }
        whiteLists2 = params.whitelist
        //if(whiteLists2 != "null"){
        //    whiteLists2 = whiteLists { it == '' ? params.whitelist : it }
        //} else{
        //    whileLists2 = "FALSE"
        //}
    }
    bambuTxDisc_out_ch = bambu(sampleIds, bamsFiles, "$params.genome", "$params.annotation", "$params.bambuPath", params.bambuParams,"$params.NDR",barcodeMaps2, whiteLists2, "$params.lowMemory")
    bambuQuant_out_ch = bambu_EM(bambuTxDisc_out_ch, "$params.genome", "$params.bambuPath", "$params.clusters", "$params.resolution")
}



