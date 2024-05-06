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
scflow_qc_options.args                 =
    "--key_colname ${params.qc_key_colname} \
    --factor_vars ${params.qc_factor_vars} \
    --min_library_size ${params.qc_min_library_size} \
    --max_library_size ${params.qc_max_library_size} \
    --min_features ${params.qc_min_features} \
    --max_features ${params.qc_max_features} \
    --max_mito ${params.qc_max_mito} \
    --min_ribo ${params.qc_min_ribo} \
    --max_ribo ${params.qc_max_ribo} \
    --min_counts ${params.qc_min_counts} \
    --min_cells ${params.qc_min_cells} \
    --drop_unmapped ${params.qc_drop_unmapped} \
    --drop_mito ${params.qc_drop_mito} \
    --drop_ribo ${params.qc_drop_ribo} \
    --nmads ${params.qc_nmads} \
    --find_singlets ${params.mult_find_singlets} \
    --singlets_method ${params.mult_singlets_method} \
    --vars_to_regress_out ${params.mult_vars_to_regress_out} \
    --pca_dims ${params.mult_pca_dims} \
    --var_features ${params.mult_var_features} \
    --doublet_rate ${params.mult_doublet_rate} \
    --dpk ${params.mult_dpk} \
    --pK ${params.mult_pK} \
    --find_cells ${params.amb_find_cells} \
    --lower ${params.amb_lower} \
    --retain ${params.amb_retain} \
    --alpha_cutoff ${params.amb_alpha_cutoff} \
    --niters ${params.amb_niters}  \
    --expect_cells ${params.amb_expect_cells} \
    --species ${params.species} "

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
    export MC_CORES=${task.cpus}

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
