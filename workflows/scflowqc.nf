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

def scflow_qc_options                  = modules['scflow_qc']

include { SCFLOW_CHECKINPUTS         } from '../modules/local/process/checkinputs'       addParams( options: scflow_checkinputs_options       )
include { SCFLOW_QC                  } from '../modules/local/process/qc'                addParams( options: scflow_qc_options                )
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process unzip_mat {
    tag "${key}"

    input:
    tuple val(key), path(mat_path)

    output:
    tuple val(key), path(mat_folder)


    script:
    """
    if [[ -d ${mat_path} ]]; then
        echo "${mat_path} is a directory"
        mv ${mat_path} mat_folder
    elif [[ -f ${mat_path} ]]; then
        echo "${mat_path} is a file"
        mkdir mat_folder && unzip ${mat_path} -d ./mat_folder
    else
        echo "${mat_path} is not valid"
        mv ${mat_path} mat_folder
        exit 1
    fi
    """

}

workflow SCFLOWQC {

    main:
    SCFLOW_CHECKINPUTS (
        ch_manifest,
        ch_input
    )

    unzip_mat (
        SCFLOW_CHECKINPUTS.out.checked_manifest.splitCsv(
            header:['key', 'filepath'],
            skip: 1, sep: '\t'
            )
        .map { row -> tuple(row.key, row.filepath) }
    )

    SCFLOW_QC (
        unzip_mat.out,
        ch_input,
        ch_ensembl_mappings
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
