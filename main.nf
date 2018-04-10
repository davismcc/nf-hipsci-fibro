#!/usr/bin/env nextflow
/*
========================================================================================
                         davismcc/nf-hipsci-fibro
========================================================================================
 davismcc/nf-hipsci-fibro Analysis Pipeline. Started 2018-04-10.
 #### Homepage / Documentation
 https://github.com/davismcc/nf-hipsci-fibro
 #### Authors
 Davis McCarthy davismcc <davis@ebi.ac.uk> - https://github.com/davismcc>
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info"""
    =========================================
     davismcc/nf-hipsci-fibro v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run davismcc/nf-hipsci-fibro --reads '*_R{1,2}.fastq.gz' -profile docker
    nextflow run davismcc/nf-hipsci-fibro -w '/hps/nobackup/hipsci/scratch/singlecell_fibroblast/nf-work' --reads '/hps/nobackup/hipsci/scratch/singlecell_fibroblast/Data/SS2_2017/22*/fastq/*_{1,2}_val_{1,2}.fq.gz' --fasta '/hps/nobackup/stegle/datasets/references/human/STAR_GRCh37.75_ERCC/GRCh37.p13.genome.ERCC92.fa' -N 'davis@ebi.ac.uk' --bams '/hps/nobackup/hipsci/scratch/singlecell_fibroblast/Data/SS2_2017/*/star/*/*realigned.bqsr.bam' -profile lsf -resume
 

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --bams                        Path to processed BAM files for variant calling
      --genome                      Name of iGenomes reference
      -profile                      Hardware config to use. lsf / docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
// params.bams = '/hps/nobackup/hipsci/scratch/singlecell_fibroblast/Data/SS2_2017/*/star/*/*.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam'

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_fastqc }

/*
 * Create a channel for processed bam files
 */
Channel
    .fromPath( params.bams )
    .ifEmpty { exit 1, "Cannot find any bams matching: ${params.bams}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
    .into { read_files_bams }



// Header log info
log.info "========================================="
log.info " davismcc/nf-hipsci-fibro v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect()
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}


/*
 * STEP 3 - Call variants: bcftools mpileup
 */

process bcftools_mpileup {
    tag "$name"
    publishDir "${params.outdir}/mpileup", mode: 'copy',
        saveAs: {"$filename"}

    input:
    set val(name), file(bam) from read_files_bams
    file fasta from fasta

    output:
    file "*vcf.gz" into mpileup_results

    script:
    """
    bcftools mpileup -E -Oz -R ${params.sites} -f ${fasta} -o "*vcf.gz" ${bam}
    """
}


/*
 * STEP 4 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


