include { FQTK } from './modules/local/fqtk/main'

params.outdir = './results'
params.readstructure_samplesheet = false
params.barcodes_samplesheet = false

process CREATE_SAMPLESHEET {
    input:
    val(demultiplexed_fastqs) // Flat list of fastq file paths
    val(baseDir)

    output:
    path('samplesheet.csv')

    script:
    /*
    Desired samplesheet.csv format:

    sample,fastq_1,fastq_2,strandedness
    A1,A1.R1.fq.gz,A1.R2.fq.gz,auto  // Paired-end
    B1,B1.R1.fq.gz,,auto            // Single-end (only R1 found)

    Pairs files based on base name before .R1 or .R2
    If baseDir is provided, paths will be relative to it.
    */
    def sampleMap = [:] // Map to store [sample_name: [R1: path, R2: path]]

    // Iterate over the flat list of fastq files
    demultiplexed_fastqs.each { fastq_path ->
        def baseName = fastq_path.baseName
        def sample_name = ""
        def read_num = ""

        // Determine sample name and read number (R1/R2)
        if (baseName.contains('.R1.')) {
            sample_name = baseName.replaceAll(/\.R1\..*/, "")
            read_num = 'R1'
        } else if (baseName.contains('.R2.')) {
            sample_name = baseName.replaceAll(/\.R2\..*/, "")
            read_num = 'R2'
        } else {
            // Handle files that don't match R1/R2 pattern if necessary
            // For now, we assume they might be single-end without the R1 marker
            // Or skip them if they are unexpected files
            // Let's assume for now it's an SE read without R1 marker
            sample_name = baseName.replaceAll(/\..*/, "") // Basic name extraction
            read_num = 'R1' // Treat as R1
        }

        if (sample_name) {
            if (!sampleMap.containsKey(sample_name)) {
                sampleMap[sample_name] = [:]
            }
            sampleMap[sample_name][read_num] = fastq_path
        }
    }

    // Build the CSV content string
    def csv_content = "sample,fastq_1,fastq_2,strandedness\n"
    sampleMap.each { sample_name, reads ->
        def r1_path_obj = reads.containsKey('R1') ? reads['R1'] : null
        def r2_path_obj = reads.containsKey('R2') ? reads['R2'] : null

        if (r1_path_obj) { // Only output if R1 exists
            def final_r1_path
            def final_r2_path

            if (baseDir && baseDir.trim()) { // Check if baseDir is a non-empty string
                // Construct the final paths using baseDir and the file name
                final_r1_path = "${baseDir}/${r1_path_obj.name}"
                final_r2_path = r2_path_obj ? "${baseDir}/${r2_path_obj.name}" : ''
            } else {
                // Use the original paths (from work directory)
                final_r1_path = r1_path_obj.toString()
                final_r2_path = r2_path_obj ? r2_path_obj.toString() : ''
            }

             csv_content += "${sample_name},${final_r1_path},${final_r2_path},auto\n"
        }
        // Optionally add logging here for cases where only R2 is found, which might indicate an issue
        // else if (r2_path) { log.warn "Found R2 without R1 for sample ${sample_name}: ${r2_path}" }
    }

    """
    cat <<EOF > samplesheet.csv
${csv_content}EOF
    """
}

workflow {

    main:
    // Read the read structure samplesheet file directly
    def readStructureFile = file(params.readstructure_samplesheet)
    if (!readStructureFile.exists()) {
        error "Read structure samplesheet not found: ${params.readstructure_samplesheet}"
    }

    // Parse the read structure samplesheet to create the list:
    // [[<fastq name: string>, <read structure: string>, <path to fastqs: path>], [example_R1.fastq.gz, 150T, ./work/98/30bc..78y/fastqs/]]
    List fqtk_readstructure_list = []
    readStructureFile.splitCsv(header: true, sep: '\t', strip: true).each { row ->
        def fastqName = file(row.filename).name
        def readStructure = row.read_structure
        def fastqPath = file(row.filename).parent
        fqtk_readstructure_list.add([fastqName, readStructure, fastqPath])
    }

    // Get the path to the barcodes samplesheet
    def barcodesSamplesheetPath = file(params.barcodes_samplesheet)
    if (!barcodesSamplesheetPath.exists()) {
        error "Barcodes samplesheet not found: ${params.barcodes_samplesheet}"
    }

    // Define metadata
    // Using the simple name of the read structure file as a basic ID
    def meta = [ id: readStructureFile.simpleName ]

    // Create the input channel for FQTK with the required structure
    // tuple val(meta), path(sample_sheet), val(fastq_readstructure_pairs)
    ch_fqtk_input = Channel.of([ meta, barcodesSamplesheetPath, fqtk_readstructure_list ])

    ch_fqtk_input.view()

    FQTK(ch_fqtk_input)

    // Collect the FASTQ pairs output by FQTK
    // Use map to discard the meta part and keep only the path part
    ch_collected_fastqs = FQTK.out.sample_fastq
        .map { _meta, fastq -> fastq }
        .collect()

    // We change the fastq files that have './work' directory paths
    // to have paths in the output (publishDir) directory instead
    // We use getAbsoluteFile and canonicalPath to do something like 'realpath'
    demuxed_basedir = 
        (new java.io.File("${params.outdir}/demultiplexed/output")).getAbsoluteFile().canonicalPath

    CREATE_SAMPLESHEET(
        ch_collected_fastqs, 
        demuxed_basedir)
}