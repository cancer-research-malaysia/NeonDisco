#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Function to read manifest file (same as yours but with debugging)
def createInputChannelFromManifestDebug(manifestPath) {
    log.info "[DEBUG] Reading manifest from: ${manifestPath}"
    
    def manifestFile = file(manifestPath)
    if (!manifestFile.exists()) {
        log.error "Manifest file not found at path: ${manifestPath}"
        exit 1
    }
    
    log.info "[DEBUG] Manifest file exists, reading contents..."

    def inputCh = channel
        .fromPath(manifestPath)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            log.info "[DEBUG] Processing row: ${row}"
            
            def sampleName = row.sampleName
            def read1 = row.read1Path
            def read2 = row.read2Path
            def sampleType = row.sampleType
            
            log.info "[DEBUG] Extracted - Sample: ${sampleName}, Type: ${sampleType}, Read1: ${read1}"
            
            // Validate required fields
            if (!sampleName || !read1 || !read2 || !sampleType) {
                log.error "[DEBUG] Invalid manifest row: ${row}. Missing fields."
                return null
            }
            
            log.info "[DEBUG] Row valid, returning tuple"
            return tuple(sampleName, [read1, read2], sampleType)
        }
        .filter { rowData -> 
            def isValid = rowData != null
            log.info "[DEBUG] Filter result: ${isValid ? 'PASS' : 'FAIL'}"
            return isValid
        }

    return inputCh
}

// Function to branch with debugging
def branchInputChannelBySampleTypeDebug(inputCh) {
    log.info "[DEBUG] Starting to branch input channel..."
    
    def branchedCh = inputCh
        .map { sample ->
            log.info "[DEBUG] Pre-branch: Sample=${sample[0]}, Type=${sample[2]}"
            return sample
        }
        .branch { sample ->
            tumor: sample[2] == 'Tumor'
            normal: sample[2] == 'Normal'
        }
    
    // Convert back to the original format
    def tumorCh = branchedCh.tumor.map { sampleName, reads, sampleType -> 
        log.info "[DEBUG] Tumor channel - emitting: ${sampleName}"
        return tuple(sampleName, reads) 
    }
    
    def normalCh = branchedCh.normal.map { sampleName, reads, sampleType -> 
        log.info "[DEBUG] Normal channel - emitting: ${sampleName}"
        return tuple(sampleName, reads) 
    }
    
    return [tumorCh, normalCh]
}

workflow {
    // Test with your manifest
    params.manifestPath = 'manifests/edgren-et-al-manifest.tsv'
    
    log.info "=== STARTING DEBUG WORKFLOW ==="
    log.info "Manifest path: ${params.manifestPath}"
    
    // First, let's just see what's in the file
    log.info "[DEBUG] Reading manifest file directly..."
    channel.fromPath(params.manifestPath)
        .splitCsv(header: true, sep: '\t')
        .subscribe { row ->
            log.info "[DEBUG RAW ROW] ${row}"
        }
    
    // Now run through the full process
    log.info "[DEBUG] Creating input channel..."
    def inputCh = createInputChannelFromManifestDebug(params.manifestPath)
    
    // Count total inputs
    inputCh
        .toList()
        .subscribe { allSamples ->
            log.info "[DEBUG] Total samples from manifest: ${allSamples.size()}"
            allSamples.each { sample ->
                log.info "[DEBUG] Sample details: Name=${sample[0]}, Type=${sample[2]}"
            }
        }
}

workflow TEST_BRANCH {
    log.info "=== TESTING BRANCH LOGIC ==="
    
    def inputCh = createInputChannelFromManifestDebug(params.manifestPath)
    
    def (tumorCh, normalCh) = branchInputChannelBySampleTypeDebug(inputCh)
    
    // Count tumor samples
    tumorCh
        .toList()
        .subscribe { tumorList ->
            log.info "[DEBUG] Total TUMOR samples after branch: ${tumorList.size()}"
            tumorList.each { sample ->
                log.info "[DEBUG TUMOR] ${sample[0]}"
            }
        }
    
    // Count normal samples
    normalCh
        .toList()
        .subscribe { normalList ->
            log.info "[DEBUG] Total NORMAL samples after branch: ${normalList.size()}"
            normalList.each { sample ->
                log.info "[DEBUG NORMAL] ${sample[0]}"
            }
        }
}

workflow TEST_MULTIMAP {
    log.info "=== TESTING MULTIMAP LOGIC ==="
    
    def inputCh = createInputChannelFromManifestDebug(params.manifestPath)
    def (tumorCh, normalCh) = branchInputChannelBySampleTypeDebug(inputCh)
    
    log.info "[DEBUG] Setting up multiMap..."
    
    def tumorSplitCh = tumorCh.multiMap { sample ->
        log.info "[DEBUG MULTIMAP] Processing: ${sample[0]}"
        forCount: sample
        forPipeline: sample
    }
    
    log.info "[DEBUG] Counting samples..."
    
    tumorSplitCh.forCount
        .toList()
        .subscribe { sampleList -> 
            log.info "[DEBUG] Final tumor sample count: ${sampleList.size()}"
            sampleList.each { sample ->
                log.info "[DEBUG COUNT] ${sample[0]}"
            }
        }
    
    tumorSplitCh.forPipeline
        .toList()
        .subscribe { sampleList ->
            log.info "[DEBUG] Samples going to pipeline: ${sampleList.size()}"
        }
}
