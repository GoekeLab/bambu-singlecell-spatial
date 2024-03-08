#! /usr/bin/env nextflow

nextflow.enable.dsl=2

params.reads = "input/sampled.fastq.gz"
params.annotation = "input/Homo_sapiens.GRCh38.91.sorted.gtf"
params.genome = "input/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
params.chemistry = "3prime-v3"
params.technology = "ONT"

params.ncore = 4

params.outdir = "output" 

process flexiplex{
	
	publishDir "$params.outdir", mode: 'copy' 
	cpus params.ncore
	maxForks params.ncore
	
	input: 
	path(fastq) 
	path(whitelist)
	
	output: 
	path ('new_reads.fastq') 
	
	script: 
	""" 

	
	if [[ $params.chemistry == 10x3v2 ]]; then 
		leftFlanking="CTACACGACGCTCTTCCGATCT"
		rightFlanking="TTTTTTTTT"
		barcodePlaceholder="????????????????"
		umiPlaceholder="??????????"
	elif [ $params.chemistry == 10x3v3 ]]; then 
		leftFlanking="CTACACGACGCTCTTCCGATCT"
		rightFlanking="TTTTTTTTT"
		barcodePlaceholder="????????????????"
		umiPlaceholder="????????????"

	elif [[ $params.chemistry = 10x5v2 ]]; then 
		leftFlanking="CTACACGACGCTCTTCCGATCT"
		rightFlanking="TTTCTTATATGGG"
		barcodePlaceholder="????????????????"
		umiPlaceholder="??????????"
		
	fi

	gunzip -c $whitelist > whitelist.txt 
	gunzip -c $fastq > reads.fastq 
	flexiplex -p $params.ncore -x \$leftFlanking -b \$barcodePlaceholder -u \$umiPlaceholder -x \$rightFlanking -f 0 reads.fastq
	python $projectDir/filter-barcodes.py --whitelist whitelist.txt --outfile my_filtered_barcode_list.txt flexiplex_barcodes_counts.txt
	flexiplex -p $params.ncore -k my_filtered_barcode_list.txt -x \$leftFlanking -b \$barcodePlaceholder -u \$umiPlaceholder -x \$rightFlanking -f 8 -e 1 reads.fastq > new_reads.fastq
	
	""" 

}

process minimap{ 

	cpus params.ncore
	maxForks params.ncore
	
	input: 
	path(newfastq)
	path(genome)

	output: 
	path ('demultiplexed.bam') 

	script:
	""" 
	if [[ $params.technology == PacBio ]]; then 
		minimap2 -ax splice:hq -t $params.ncore -d ref.mmi $genome
		minimap2 -ax splice:hq -t $params.ncore -a ref.mmi $newfastq > demultiplexed.sam  
	
	elif [[ $params.technology == ONT ]]; then 
		minimap2 -ax splice -k14 -t $params.ncore -d ref.mmi $genome
		minimap2 -ax splice -k14 -t $params.ncore -a ref.mmi $newfastq > demultiplexed.sam  

	fi

	samtools sort -@ $params.ncore demultiplexed.sam -o demultiplexed.bam 
	samtools index -@ $params.ncore demultiplexed.bam 

	rm demultiplexed.sam
	rm ref.mmi 
	"""
}

process bambu_discovery{ 
    publishDir "$PWD/$params.outdir", mode: 'copy' 

	cpus params.ncore 
	maxForks params.ncore
	
	input: 
	path(bam)
	path(genome)
	path(annotation)

	output: 
	path ('readClassFile.rds')
	path ('extended_annotations.rds') 
    	path ('extended_annotations_NDR1.rds') 
    	path ('quantData.rds')

	script:
	""" 
	#!/usr/bin/env Rscript
	library(devtools)
	load_all("$PWD/bambu")
    	library(bambu)

	annotations <- prepareAnnotations("$annotation")

	# Transcript discovery and generate readGrgList for each cell
    readClassFile = bambu(reads = "$bam", annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = TRUE, assignDist = FALSE, lowMemory = TRUE, yieldSize = 10000000)
    print("readClassDone")
    saveRDS(readClassFile, "readClassFile.rds")
    extendedAnno = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    saveRDS(extendedAnno, "extended_annotations.rds")
    extendedAnno_NDR1 = bambu(reads = readClassFile, annotations = annotations, genome = "$genome", NDR = 1, ncore = $params.ncore, discovery = TRUE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, assignDist = FALSE)
    saveRDS(extendedAnno_NDR1, "extended_annotations_NDR1.rds")
    print("extendAnno done")
    rm(extendedAnno_NDR1)
    rm(annotations)
    se = bambu(reads = readClassFile, annotations = extendedAnno, genome = "$genome", ncore = $params.ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE)
    saveRDS(se, "quantData.rds")
	"""
}

process bambu_quant{

	publishDir "$PWD/$params.outdir", mode: 'copy'

	cpus params.ncore
        maxForks params.ncore

        input:
        path(readClassFile)
	path(extendedAnno)
	path(extendedAnno_NDR1)
	path(quantData)
        path(genome)

        output:
        path ('se.rds')

        script:
	"""
        #!/usr/bin/env Rscript
        library(devtools)
        load_all("$PWD/bambu")
        library(bambu)
	extendedAnno <- readRDS("$extendedAnno")
	quantData = readRDS("$quantData")
        se = bambu(reads = "test.rds", annotations = extendedAnno, genome = "$genome", quantData = quantData, assignDist = FALSE, ncore = $params.ncore, discovery = FALSE, quant = TRUE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE))

        saveRDS(se, "se.rds")
        writeBambuOutput(se, path = ".", prefix = "sparse_")
        """

}

// This is the workflow to execute the process 
workflow {
	flexiplex_out_ch = flexiplex("$PWD/$params.reads", "$PWD/$params.whitelist")
	minimap_out_ch = minimap(flexiplex_out_ch, "$PWD/$params.genome")
	bambuTxDisc_out_ch = bambu_discovery(minimap_out_ch, "$PWD/$params.genome", "$PWD/$params.annotation")
	bambuQuant_out_ch = bambu_quant(bambuTxDisc_out_ch, "$PWD/$params.genome")
}

