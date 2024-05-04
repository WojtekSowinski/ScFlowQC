/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.manifest) { ch_manifest = file(params.manifest, checkIfExists: true) }
if (params.input) { ch_input = file(params.input, checkIfExists: true)}
if (params.ensembl_mappings) { ch_ensembl_mappings = file(params.ensembl_mappings, checkIfExists: false) }


def modules = params.modules.clone()

def scflow_checkinputs_options = modules['scflow_checkinputs']
scflow_checkinputs_options.args = ''

def scflow_qc_options = modules['scflow_qc']
scflow_qc_options.args = '' // TODO: define QC options???

include { SCFLOW_CHECKINPUTS         } from '../modules/local/process/checkinputs'       addParams( options: scflow_checkinputs_options       )
include { SCFLOW_QC                  } from '../modules/local/process/qc'                addParams( options: scflow_qc_options                )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCFLOWQC {

    main:
    SCFLOW_CHECKINPUTS (
        ch_manifest,
        ch_input
    )

    SCFLOW_QC (
        SCFLOW_CHECKINPUTS.out.checked_manifest.splitCsv(
            header:['key', 'filepath'],
            skip: 1, sep: '\t'
            )
        .map { row -> tuple(row.key, row.filepath) },
        ch_input,
        ch_ensembl_mappings
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
