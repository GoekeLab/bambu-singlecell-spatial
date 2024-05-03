#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.chemistry = "custom" //"10x3v3" "10x3v2" "10x5v2"
params.technology = 'ONT' //"ONT" "PacBio"
params.whitelist = "nowhitelistfile"
params.bambuPath = "bambu"
params.lowMemory = false
params.ncore = 4
params.spatial = false

params.NDR = "NULL"
params.cleanReads = "TRUE"
params.keepChimericReads = false
params.deduplicateUMIs = "TRUE"
params.bambuParams = tuple(params.cleanReads, params.keepChimericReads, params.deduplicateUMIs)
params.barcodeMap = "TRUE"
//params.reads = false
//params.bams = false

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
    if [[ "$whitelist" == "nowhitelistfile" ]]; then
        python /mnt/software/filter-barcodes.py --outfile my_filtered_barcode_list.txt flexiplex_barcodes_counts.txt 
    else
        gunzip -c $whitelist > whitelist.txt 
	    python /mnt/software/filter-barcodes.py --outfile my_barcode_list.txt flexiplex_barcodes_counts.txt 
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

process bambu_discovery{ 
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

	output: 
    tuple val ('combined'), path ('*readClassFile.rds'), path ('*quantData.rds')
	path ('*extended_annotations.rds') 
    path ('*extended_annotations_NDR1.rds') 
    path ('*extended_annotations_NDR1.gtf') 

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

    if("$bambuPath" == "bambu") {
        library(bambu)
    } else {
        library(devtools)
        load_all("$bambuPath")
    }

	annotations <- prepareAnnotations("$annotation")

	# Transcript discovery and generate readGrgList for each cell
    readClassFile = bambu(reads = samples, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = barcode_maps, verbose = TRUE, assignDist = FALSE, lowMemory = TRUE, yieldSize = 10000000, sampleNames = ids, cleanReads = as.logical($cleanReads), dedupUMI = as.logical($deduplicateUMIs))
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
    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE)
    saveRDS(se, paste0(runName, "_quantData.rds"))
    #write(runName, "runName.txt")
	"""
}

process bambu_readClassConstruction{
    cpus params.ncore
    maxForks 1

    input:
	tuple val(runName),path(bam)
	path(genome)
	path(annotation)
    val(bambuPath)
    tuple val(cleanReads),val(keepChimericReads), val(deduplicateUMIs)

	output: 
    tuple val(runName), path ('*readClassFile.rds')

    script:
	""" 
	#!/usr/bin/env Rscript

    if("$bambuPath" == "bambu") {
        library(bambu)
    } else {
        library(devtools)
        load_all("$bambuPath")
    }
    annotations <- prepareAnnotations("$annotation")
    readClassFile = bambu(reads = "$bam", annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = TRUE, assignDist = FALSE, lowMemory = TRUE, yieldSize = 10000000, cleanReads = as.logical($cleanReads), dedupUMI = as.logical($deduplicateUMIs))
    saveRDS(readClassFile, paste0("$runName", "_readClassFile.rds"))
    """

}

process bambu_extend{
    publishDir "$params.outdir", mode: 'copy'
    
    cpus params.ncore
    maxForks params.ncore

    input:
    val(readClassFile)
    path(genome)
    path(annotation)
    val(bambuPath)
    val(NDR)

	output: 
    path ('*extended_annotations.rds') 
    path ('*extended_annotations_NDR1.rds') 
    path ('extended_annotations_NDR1.gtf')

    script:
	""" 
	#!/usr/bin/env Rscript

    samples = "$readClassFile"
    samples = gsub("[][]","", gsub(' ','', samples))
    samples = unlist(strsplit(samples, ','))
    print(samples)

    if("$bambuPath" == "bambu") {
        library(bambu)
    } else {
        library(devtools)
        load_all("$bambuPath")
    }
    annotations <- prepareAnnotations("$annotation")
    if(isFALSE($NDR)){
        extendedAnno = bambu(reads = samples, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    } else{
        extendedAnno = bambu(reads = samples, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE, NDR = $NDR)
    }
    saveRDS(extendedAnno, "_extended_annotations.rds")
    extendedAnno_NDR1 = bambu(reads = samples, annotations = annotations, genome = "$genome", NDR = 1, ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    saveRDS(extendedAnno_NDR1, "_extended_annotations_NDR1.rds")
    writeToGTF(extendedAnno_NDR1, "extended_annotations_NDR1.gtf")
    """

}

process bambu_assignDist{
    cpus params.ncore
    maxForks params.ncore

    input:
    tuple val(runName),path(readClassFile)
    path(extendedAnno)
    path(extendedAnno_NDR1)
    path(genome)
    val(bambuPath)

	output: 
    tuple val(runName), path ('*quantData.rds'), emit: quantData

    script:
	""" 
	#!/usr/bin/env Rscript

    if("$bambuPath" == "bambu") {
        library(bambu)
    } else {
        library(devtools)
        load_all("$bambuPath")
    }
    extendedAnno <- readRDS("$extendedAnno")
    se = bambu(reads = "$readClassFile", annotations = extendedAnno, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE)
    saveRDS(se, paste0("$runName", "_quantData.rds"))
    """
}

process bambu_quant{

	publishDir "$params.outdir", mode: 'copy'

	cpus params.ncore
    maxForks 1

    input:
    tuple val(sample),path(readClassFile), path(quantData)
    path(extendedAnno)
    path(extendedAnno_NDR1)
    path(extended_annotations_NDR1_gtf) 
    path(genome)
    val(bambuPath)
    val(whitelist)

    output:
    path ('*.rds')
    path ('*.mtx')
    path ('*.gtf')
    path ('*.tsv')

    script:
    """
    #!/usr/bin/env Rscript
    if("$bambuPath" == "bambu") {
        library(bambu)
    } else {
        library(devtools)
        load_all("$bambuPath")
    }
    if(".txt" %in% "$sample"){runName = readLines("$sample")
    } else{runName = "$sample"}
    
    extendedAnno <- readRDS("$extendedAnno")
    quantData = readRDS("$quantData")
    if("$whitelist" == FALSE){
        se = bambu(reads = "test.rds", annotations = extendedAnno, genome = "$genome", quantData = quantData, assignDist = FALSE, ncore = $params.ncore, discovery = FALSE, quant = TRUE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE))
    } else{
        se = bambu(reads = "test.rds", annotations = extendedAnno, genome = "$genome", quantData = quantData, assignDist = FALSE, ncore = $params.ncore, spatial = "$whitelist", discovery = FALSE, quant = TRUE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE))
    }
    saveRDS(se, paste0(runName, "_se.rds"))

    writeBambuOutput(se, path = ".", prefix = paste0(runName, "_sparse_"))
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
                    | map { row-> tuple(row.sample, file(row.fastq), row.chemistry, row.technology, row.whitelist) }
            flexiplex_out_ch = flexiplex(readsChannel)
            minimap_out_ch = minimap(flexiplex_out_ch, "$params.genome")
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
            //params.bams = minimap_out_ch
        }

    }

    if(!params.reads) {
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
        if(params.lowMemory){
            //bambuTxDisc_out_ch = bambu_discovery(bamsChannel, "$params.genome", "$params.annotation")
            //bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$params.genome")
            readClassConstruction_out_ch = bambu_readClassConstruction(bamsChannel, "$params.genome", "$params.annotation", "$params.bambuPath", params.bambuParams)
            extend_out_ch = bambu_extend(readClassConstruction_out_ch.collect{it[1]},"$params.genome","$params.annotation", "$params.bambuPath", "$params.NDR")
            assignDist_out_ch = bambu_assignDist(readClassConstruction_out_ch, extend_out_ch, "$params.genome", "$params.bambuPath")
            readClassConstruction_out_ch
                .join(assignDist_out_ch, by:0)
                .set { bambuQuant_inputs }
            bambuQuant_out_ch =  bambu_quant(bambuQuant_inputs, extend_out_ch, "$params.genome", "$params.bambuPath", "$params.spatial")
        }
        else{
            sampleIds = bamsChannel.collect{it[0]}
            bamsFiles = bamsChannel.collect{it[1]}
            barcodeMaps = bamsChannel.collect{it[2]}
            whiteLists = bamsChannel.collect{it[3]}

            barcodeMaps2 = barcodeMaps.map { it == null ? it : "TRUE" }

            if(params.whitelist == "nowhitelistfile" & whiteLists.any { it == '' }){
                echo "Missing whitelist entries in samplesheet or no --whitelist provided"
                exit 1
            }
            whiteLists2 = params.whitelist
            if(whiteLists2 != "nowhitelistfile"){
                whiteLists2 = whiteLists { it == '' ? params.whitelist : it }
            }


            bambuTxDisc_out_ch = bambu_discovery(sampleIds, bamsFiles, "$params.genome", "$params.annotation", "$params.bambuPath", params.bambuParams,"$params.NDR",barcodeMaps2)
            if(params.spatial){
                bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$params.genome", "$params.bambuPath", whiteLists2)
            }
            else{
                bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$params.genome", "$params.bambuPath", "FALSE")
            }
        }

    }
    else{ //When preprocessing steps were run
        if(params.lowMemory){
            //bambuTxDisc_out_ch = bambu_discovery(params.bams, "$params.genome", "$params.annotation")
            //bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$params.genome")
            readClassConstruction_out_ch = bambu_readClassConstruction(params.bams, "$params.genome", "$params.annotation", "$params.bambuPath", params.bambuParams)
            extend_out_ch = bambu_extend(readClassConstruction_out_ch.collect{it[1]},"$params.genome","$params.annotation", "$params.bambuPath", "$params.NDR")
            assignDist_out_ch = bambu_assignDist(readClassConstruction_out_ch, extend_out_ch, "$params.genome", "$params.bambuPath")
            readClassConstruction_out_ch
                .join(assignDist_out_ch, by:0)
                .set { bambuQuant_inputs }
            bambuQuant_out_ch =  bambu_quant(bambuQuant_inputs, extend_out_ch, "$params.genome", "$params.bambuPath", "$params.spatial")        }
        else{
            sampleIds = minimap_out_ch.collect{it[0]}
            bamsFiles = minimap_out_ch.collect{it[1]}
            bambuTxDisc_out_ch = bambu_discovery(sampleIds, bamsFiles, "$params.genome", "$params.annotation", "$params.bambuPath", params.bambuParams,"$params.NDR", params.barcodeMap)
            if(params.spatial){
                bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$params.genome", "$params.bambuPath", params.whitelist)
            }
            else{
                bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$params.genome", "$params.bambuPath", "FALSE")
            }
        }
    }
}



