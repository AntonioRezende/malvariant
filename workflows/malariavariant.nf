/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMalariavariant.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                         } from '../modules/nf-core/fastqc/main'
include { MULTIQC                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS    } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FASTP                          } from '../modules/nf-core/fastp/main'
include { PICARD_FASTQTOSAM              } from '../modules/nf-core/picard/fastqtosam/main'
include { BWA_INDEX                      } from '../modules/nf-core/bwa/index/main' 
include { BWA_MEM                        } from '../modules/nf-core/bwa/mem/main'
include { GATK4_MERGEBAMALIGNMENT        } from '../modules/nf-core/gatk4/mergebamalignment/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { SAMTOOLS_FAIDX                 } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_FILTERMUTECTCALLS        } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { GATK4_SELECTVARIANTS           } from '../modules/nf-core/gatk4/selectvariants/main'
//include { GATK4_VARIANTFILTRATION        } from '../modules/nf-core/gatk4/variantfiltration/main'
include { BCFTOOLS_VIEW                  } from '../modules/nf-core/bcftools/view/main'
//include { SNPEFF_SNPEFF                  } from '../modules/nf-core/snpeff/snpeff/main'
include { GATK4_GENOMICSDBIMPORT         } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS            } from '../modules/nf-core/gatk4/genotypegvcfs/main'

include {PICARD_MARKILLUMINAADAPTERS     } from '../modules/local/markadapters'
include {PICARD_SAMTOFASTQ               } from '../modules/local/samtofastq'
include {GATK4_MUTECT2                   } from '../modules/local/mutect2'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MALARIAVARIANT {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Fastp trimming
    //
    FASTP (
        INPUT_CHECK.out.reads,
        [],
        false,
        false,
        false
    )
    //FASTP.out.reads.view()
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    PICARD_FASTQTOSAM(
        FASTP.out.reads
    )
    //PICARD_FASTQTOSAM.out.bam.view()
    ch_versions = ch_versions.mix(PICARD_FASTQTOSAM.out.versions.first())
    

    PICARD_MARKILLUMINAADAPTERS{
        PICARD_FASTQTOSAM.out.bam
    }
    //PICARD_MARKILLUMINAADAPTERS.out.bam.view()
    ch_versions = ch_versions.mix(PICARD_MARKILLUMINAADAPTERS.out.versions.first())

    PICARD_SAMTOFASTQ{
        PICARD_MARKILLUMINAADAPTERS.out.bam
    }
    //PICARD_SAMTOFASTQ.out.reads.view()
    ch_versions = ch_versions.mix(PICARD_SAMTOFASTQ.out.versions.first())

    reference = [
        [ id:'reference', single_end:true ],
        file(params.fasta)
    ]

    bedregions =[
         [ id:'reference', single_end:true ],
        file(params.bedfile)
    ]

    //cachesnpEff =[
    //    [id:'reference', single_end:true ],
    //    file("/ARezende/singularityNF")
    //]
    //cachesnpEff.view()
    
    BWA_INDEX{
        reference
    }
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    GATK4_CREATESEQUENCEDICTIONARY{
        reference
    }
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions.first())
    
    BWA_MEM(
        PICARD_SAMTOFASTQ.out.reads,
        BWA_INDEX.out.index,
        reference,
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    BAMS = BWA_MEM.out.bam.join(PICARD_FASTQTOSAM.out.bam)
    //BAMS.view()

    GATK4_MERGEBAMALIGNMENT(
        BAMS,
        reference,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    ch_versions = ch_versions.mix(GATK4_MERGEBAMALIGNMENT.out.versions.first())


    SAMTOOLS_FAIDX(
        reference
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    GATK4_MUTECT2(
        GATK4_MERGEBAMALIGNMENT.out.bam,
        reference,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    //GATK4_MUTECT2.out.vcf.view()
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    tempCH = GATK4_MUTECT2.out.vcf.join(GATK4_MUTECT2.out.tbi)
    nofilteredVCF = tempCH.join(GATK4_MUTECT2.out.stats)

    GATK4_FILTERMUTECTCALLS(
        nofilteredVCF,
        reference,
        SAMTOOLS_FAIDX.out.fai,
        GATK4_CREATESEQUENCEDICTIONARY.out.dict
    )
    //GATK4_FILTERMUTECTCALLS.out.vcf.view()
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())


    filteredVCF = GATK4_FILTERMUTECTCALLS.out.vcf.join(GATK4_FILTERMUTECTCALLS.out.tbi)

    GATK4_SELECTVARIANTS(
        filteredVCF
    )
    //GATK4_SELECTVARIANTS.out.vcf
    ch_versions = ch_versions.mix(GATK4_SELECTVARIANTS.out.versions.first())

    filteredSNPs = GATK4_SELECTVARIANTS.out.vcf.join(GATK4_SELECTVARIANTS.out.tbi)

    //GATK4_VARIANTFILTRATION(
    //    filteredSNPs,
    //    reference,
    //   SAMTOOLS_FAIDX.out.fai,
    //    GATK4_CREATESEQUENCEDICTIONARY.out.dict
    //)

    BCFTOOLS_VIEW(
        filteredSNPs,
        bedregions
    )
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions.first())
    //bcfFiltered = BCFTOOLS_VIEW.out.vcf.join(BCFTOOLS_VIEW.out.tbi)
    //vcfsfrombcf=bcfFiltered.map{it[1]}.collect().toList()//.view()
    //tbifrombcf=bcfFiltered.map{it[2]}.collect().toList()//.view()
    
    //allsamplesVCF= vcfsfrombcf.merge(tbifrombcf)
    //.map{vcf,tbi -> 
    //    [[id: "VCF"], vcf, tbi]
    //}

    //allsamplesVCF.view()

    //SNPEFF_SNPEFF(
    //    BCFTOOLS_VIEW.out.vcf,
    //    "Plasmodium_falciparum",
    //    "/ARezende/singularityNF"
    //)

    //GATK4_GENOMICSDBIMPORT(
    //    allsamplesVCF,
    //    bedregions,
    //    [],
    //    [],
    //    []
    //)

    //GATK4_GENOTYPEGVCFS(
    //    GATK4_GENOMICSDBIMPORT.out.genomicsdb,
    //    reference,
    //    SAMTOOLS_FAIDX.out.fai,
    //    GATK4_CREATESEQUENCEDICTIONARY.out.dict
    //)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMalariavariant.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMalariavariant.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKILLUMINAADAPTERS.out.stats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(GATK4_FILTERMUTECTCALLS.out.stats.collect{it[1]}.ifEmpty([]))
    //ch_multiqc_files = ch_multiqc_files.mix(PICARD_FASTQTOSAM.out.json.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
