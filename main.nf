#!/usr/bin/env nextflow

import groovy.json.*

params.input = false
params.fs = false
params.bids = false
params.bids_config = false
params.help = false
params.dti_shells = false
params.fodf_shells = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["clean_bids":"$params.clean_bids",
                "sh_fitting":"$params.sh_fitting",
                "sh_fitting_basis":"$params.sh_fitting_basis",
                "sh_fitting_order":"$params.sh_fitting_order",
                "b0_thr_extract_b0":"$params.b0_thr_extract_b0",
                "dwi_shell_tolerance":"$params.dwi_shell_tolerance",
                "dilate_b0_mask_prelim_brain_extraction":"$params.dilate_b0_mask_prelim_brain_extraction",
                "bet_prelim_f":"$params.bet_prelim_f",
                "run_dwi_denoising":"$params.run_dwi_denoising",
                "extent":"$params.extent",
                "run_gibbs_correction": "$params.run_gibbs_correction",
                "run_topup":"$params.run_topup",
                "encoding_direction":"$params.encoding_direction",
                "readout":"$params.readout",
                "run_eddy":"$params.run_eddy",
                "eddy_cmd":"$params.eddy_cmd",
                "bet_topup_before_eddy_f":"$params.bet_topup_before_eddy_f",
                "use_slice_drop_correction":"$params.use_slice_drop_correction",
                "bet_dwi_final_f":"$params.bet_dwi_final_f",
                "fa_mask_threshold":"$params.fa_mask_threshold",
                "run_resample_dwi":"$params.run_resample_dwi",
                "dwi_resolution":"$params.dwi_resolution",
                "dwi_interpolation":"$params.dwi_interpolation",
                "max_dti_shell_value":"$params.max_dti_shell_value",
                "min_fodf_shell_value":"$params.min_fodf_shell_value",
                "run_t1_denoising":"$params.run_t1_denoising",
                "run_resample_t1":"$params.run_resample_t1",
                "t1_resolution":"$params.t1_resolution",
                "t1_interpolation":"$params.t1_interpolation",
                "number_of_tissues":"$params.number_of_tissues",
                "fa":"$params.fa",
                "min_fa":"$params.min_fa",
                "roi_radius":"$params.roi_radius",
                "set_frf":"$params.set_frf",
                "manual_frf":"$params.manual_frf",
                "mean_frf":"$params.mean_frf",
                "sh_order":"$params.sh_order",
                "basis":"$params.basis",
                "fodf_metrics_a_factor":"$params.fodf_metrics_a_factor",
                "relative_threshold":"$params.relative_threshold",
                "max_fa_in_ventricle":"$params.max_fa_in_ventricle",
                "min_md_in_ventricle":"$params.min_md_in_ventricle",
                "run_pft_tracking":"$params.run_pft_tracking",
                "pft_seeding_mask_type":"$params.pft_seeding_mask_type",
                "pft_fa_seeding_mask_threshold":"$params.pft_fa_seeding_mask_threshold",
                "pft_algo":"$params.pft_algo",
                "pft_seeding":"$params.pft_seeding",
                "pft_nbr_seeds":"$params.pft_nbr_seeds",
                "pft_step":"$params.pft_step",
                "pft_theta":"$params.pft_theta",
                "pft_min_len":"$params.pft_min_len",
                "pft_max_len":"$params.pft_max_len",
                "pft_compress_streamlines":"$params.pft_compress_streamlines",
                "pft_compress_value":"$params.pft_compress_value",
                "local_seeding_mask_type":"$params.local_seeding_mask_type",
                "local_fa_seeding_mask_threshold":"$params.local_fa_seeding_mask_threshold",
                "local_tracking_mask_type":"$params.local_tracking_mask_type",
                "local_fa_tracking_mask_threshold":"$params.local_fa_tracking_mask_threshold",
                "run_local_tracking":"$params.run_local_tracking",
                "local_compress_streamlines":"$params.local_compress_streamlines",
                "pft_random_seed":"$params.pft_random_seed",
                "local_algo":"$params.local_algo",
                "local_seeding":"$params.local_seeding",
                "local_nbr_seeds":"$params.local_nbr_seeds",
                "local_step":"$params.local_step",
                "local_theta":"$params.local_theta",
                "local_sfthres":"$params.local_sfthres",
                "local_sfthres_init":"$params.local_sfthres_init",
                "local_min_len":"$params.local_min_len",
                "local_max_len":"$params.local_max_len",
                "local_compress_value":"$params.local_compress_value",
                "local_random_seed":"$params.local_random_seed",
                "cpu_count":"$cpu_count",
                "template_t1":"$params.template_t1",
                "processes_brain_extraction_t1":"$params.processes_brain_extraction_t1",
                "processes_denoise_dwi":"$params.processes_denoise_dwi",
                "processes_denoise_t1":"$params.processes_denoise_t1",
                "processes_eddy":"$params.processes_eddy",
                "processes_fodf":"$params.processes_fodf",
                "processes_registration":"$params.processes_registration"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "TractoFlow pipeline"
log.info "==================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.dti_shells){
    log.info "DTI shells extracted: $params.dti_shells"
}
else{
  log.info "Max DTI shell extracted: $params.max_dti_shell_value"
}

if (params.fodf_shells){
    log.info "DTI shells extracted: $params.fodf_shells"
}
else{
  log.info "Min FODF shell extracted: $params.min_fodf_shell_value"
}


labels_for_reg = Channel.empty()
if (params.input && !(params.bids && params.bids_config)){
    log.info "Input: $params.input"
    root = file(params.input)
    data = Channel
        .fromFilePairs("$root/**/*{bval,bvec,dwi.nii.gz,t1.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}

    labels_for_reg = Channel
            .fromFilePairs("$root/**/*{aparc+aseg.nii.gz,wmparc.nii.gz}",
                           size: 2,
                           maxDepth:1,
                           flat: true) {it.parent.name}

    data
        .map{[it, params.readout, params.encoding_direction].flatten()}
        .into{in_data; check_subjects_number}

    Channel
    .fromPath("$root/**/*rev_b0.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{rev_b0_for_prepare_topup; check_simple_rev_b0}
}
else if (params.bids || params.bids_config){
    if (!params.bids_config) {
        log.info "Input BIDS: $params.bids"
        if (params.participants_label) {
            Integer start = workflow.commandLine.indexOf("participants_label") + "participants_label".length();
            participant_cleaned = workflow.commandLine.substring(start, workflow.commandLine.indexOf("--", start) == -1 ? workflow.commandLine.length() : workflow.commandLine.indexOf("--", start)).replace("=", "").replace("\'", "")
            log.info "Participants: $participant_cleaned"
        }
        if (params.fs) {
            Integer start = workflow.commandLine.indexOf("fs") + "fs".length();
            freesurfer_path = workflow.commandLine.substring(start, workflow.commandLine.indexOf("--", start) == -1 ? workflow.commandLine.length() : workflow.commandLine.indexOf("--", start)).replace("=", "").replace("\'", "")
            log.info "Freesurfer path: $freesurfer_path"
        }
        log.info "Clean_bids: $params.clean_bids"
        log.info ""

        bids = file(params.bids)

        process Read_BIDS {
            publishDir = params.Read_BIDS_Publish_Dir
            scratch = false
            stageInMode = 'symlink'
            tag = {"Read_BIDS"}
            errorStrategy = { task.attempt <= 3 ? 'retry' : 'terminate' }

            input:
            file(bids_folder) from bids

            output:
            file "tractoflow_bids_struct.json" into bids_struct

            script:
            participants_flag =\
            params.participants_label ? '--participants_label ' + participant_cleaned : ""
            fs_flag =\
            params.fs ? '--fs ' + freesurfer_path : ""

            clean_flag = params.clean_bids ? '--clean ' : ''

            """
            scil_validate_bids.py $bids_folder tractoflow_bids_struct.json\
                --readout $params.readout $participants_flag $clean_flag $fs_flag
            """
        }
    }

    else {
        log.info "BIDS config: $params.bids_config"
        config = file(params.bids_config)
        bids_struct = Channel.from(config)
    }

    ch_in_data = Channel.create()
    ch_sid_rev_dwi = Channel.create()
    ch_complex_rev_b0 = Channel.create()
    ch_simple_rev_b0 = Channel.create()
    labels_for_reg = Channel.create()

    bids_struct.map{it ->
    jsonSlurper = new JsonSlurper()
        data = jsonSlurper.parseText(it.getText())
        for (item in data){
            sid = "sub-" + item.subject

            if (item.session){
                sid += "_ses-" + item.session
            }

            if (item.run){
                sid += "_run-" + item.run
            }
            for (key in item.keySet()){
                if(item[key] == 'todo'){
                    error "Error ~ Please look at your tractoflow_bids_struct.json " +
                    "in Read_BIDS folder.\nPlease fix todo fields and give " +
                    "this file in input using --bids_config option instead of" +
                    "using --bids."
                }
                else if (item[key] == 'error_readout'){
                    error "Error ~ Please look at your tractoflow_bids_struct.json " +
                    "in Read_BIDS folder.\nPlease fix error_readout fields. "+
                    "This error indicate that readout time looks wrong.\n"+
                    "Please correct the value or remove the subject in the json and " +
                    "give the updated file in input using --bids_config option instead of" +
                    "using --bids."
                }
            }
            sub = [sid, file(item.bval), file(item.bvec), file(item.dwi), "_",
                   file(item.t1), item.TotalReadoutTime, item.DWIPhaseEncodingDir[0]]
            ch_in_data.bind(sub)

            if(item.rev_topup) {
                if(item.topup) {
                  sub_complex_rev_b0 = [sid, file(item.rev_topup), file(item.topup)]
                  ch_complex_rev_b0.bind(sub_complex_rev_b0)
                }
                else{
                  sub_simple_rev_b0 = [sid, file(item.rev_topup)]
                  ch_simple_rev_b0.bind(sub_simple_rev_b0)
                }
            }

            if(item.wmparc) {
                sub_labels_for_reg = [sid, file(item.aparc_aseg), file(item.wmparc)]
                labels_for_reg.bind(sub_labels_for_reg)
            }

            if(item.rev_dwi){
              ch_rev_in_data = [sid, file(item.rev_bval), file(item.rev_bvec), file(item.rev_dwi), "_rev_",
                                file(item.t1), item.TotalReadoutTime, item.DWIPhaseEncodingDir[0]]
              ch_sid_rev_dwi.bind([sid])
              ch_in_data.bind(ch_rev_in_data)
            }
        }
        ch_sid_rev_dwi.close()
        ch_in_data.close()
        ch_simple_rev_b0.close()
        ch_complex_rev_b0.close()
        labels_for_reg.close()
    }

    ch_sid_rev_dwi.into{sid_rev_dwi_included; sid_rev_dwi_excluded; check_rev_number}
    ch_in_data.into{in_data; check_subjects_number}

    ch_simple_rev_b0.into{rev_b0_for_prepare_topup; check_simple_rev_b0}
    ch_complex_rev_b0.into{complex_rev_b0_for_topup; check_complex_rev_b0}
}
else {
    error "Error ~ Please use --input, --bids or --bids_config for the input data."
}

check_subjects_number.map{[it[0]]}.unique().set{unique_subjects_number}

if (params.sh_fitting && !params.sh_fitting_shells){
    error "Error ~ Please set the SH fitting shell to use."
}

if (params.pft_seeding_mask_type != "wm" && params.pft_seeding_mask_type != "interface" && params.pft_seeding_mask_type != "fa"){
    error "Error ~ --pft_seeding_mask_type can only take wm, interface or fa. Please select one of these choices"
}

if (params.local_seeding_mask_type != "wm" && params.local_seeding_mask_type != "fa"){
    error "Error ~ --local_seeding_mask_type can only take wm or fa. Please select one of these choices"
}

if (params.local_tracking_mask_type != "wm" && params.local_tracking_mask_type != "fa"){
    error "Error ~ --local_tracking_mask_type can only take wm or fa. Please select one of these choices"
}

if (params.local_algo != "det" && params.local_algo != "prob"){
    error "Error ~ --local_algo can only take det or prob. Please select one of these choices"
}

if (params.pft_algo != "det" && params.pft_algo != "prob"){
    error "Error ~ --pft_algo can only take det or prob. Please select one of these choices"
}

if (params.local_seeding != "nt" && params.local_seeding != "npv"){
    error "Error ~ --local_seeding can only take nt or npv. Please select one of these choices"
}

if (params.pft_seeding != "nt" && params.pft_seeding != "npv"){
    error "Error ~ --pft_seeding can only take nt or npv. Please select one of these choices"
}

if (params.run_pft_tracking && workflow.profile.contains("ABS")){
    error "Error ~ PFT tracking cannot be run with Atlas Based Segmentation (ABS) profile"
}

if (params.bids && workflow.profile.contains("ABS") && !params.fs){
    error "Error ~ --bids parameter cannot be run with Atlas Based Segmentation (ABS) profile"
}

(dwi, gradients, t1, readout_encoding) = in_data
    .map{sid, bvals, bvecs, dwi, rev, t1, readout, encoding -> [tuple(sid, dwi, rev),
                                        tuple(sid, bvals, bvecs),
                                        tuple(sid, t1),
                                        tuple(sid, readout, encoding)]}
    .separate(4)

t1.unique()
    .into{t1_for_denoise; t1_for_test_denoise}

check_complex_rev_b0.concat(check_simple_rev_b0).count().into{rev_b0_counter; number_rev_b0_for_compare}

unique_subjects_number.count().into{ number_subj_for_null_check; number_subj_for_compare}

check_rev_number.count().set{number_rev_dwi}

if (params.eddy_cmd == "eddy_openmp"){
number_rev_dwi
    .subscribe{a -> if (a>0)
    error "Error ~ You have some subjects with a reverse encoding DWI. You MUST add -profile use_cuda with a GPU environnement to be able to analyse these data."}
}

if (!params.run_topup || !params.run_eddy){
number_rev_dwi
    .subscribe{a -> if (a>0)
    error "Error ~ You have some subjects with a reverse encoding DWI. You MUST run topup and eddy with this kind of acquisition."}
}

number_subj_for_null_check
.subscribe{a -> if (a == 0)
    error "Error ~ No subjects found. Please check the naming convention, your --input path or your BIDS folder."}

if (params.set_frf && params.mean_frf){
    error "Error ~ --set_frf and --mean_frf are activated. Please choose only one of these options. "
}

if (params.run_topup){
number_subj_for_compare
    .concat(number_rev_b0_for_compare)
    .toList()
    .subscribe{a, b -> if (a != b && b > 0)
    error "Error ~ Some subjects have a reversed phase encoded b=0 and others don't.\n" +
          "Please be sure to have the same acquisitions for all subjects."}
}

dwi.into{dwi_for_prelim_bet; dwi_for_denoise; dwi_for_test_denoise}

if (params.pft_random_seed instanceof String){
    pft_random_seed = params.pft_random_seed?.tokenize(',')
}
else{
    pft_random_seed = params.pft_random_seed
}

if (params.local_random_seed instanceof String){
    local_random_seed = params.local_random_seed?.tokenize(',')
}
else{
    local_random_seed = params.local_random_seed
}

gradients
    .into{gradients_for_prelim_bet; gradients_for_eddy;
          gradients_for_prepare_topup;
          gradients_for_prepare_dwi_for_eddy;
          gradients_for_eddy_topup; gradients_for_test_eddy_topup}

readout_encoding
    .into{readout_encoding_for_topup; readout_encoding_for_eddy;
          readout_encoding_for_eddy_topup}

dwi_for_prelim_bet
    .join(gradients_for_prelim_bet)
    .set{dwi_gradient_for_prelim_bet}

process README {
    cpus 1
    publishDir = params.Readme_Publish_Dir
    tag = "README"

    output:
    file "readme.txt"

    script:
    String list_options = new String();
    for (String item : params) {
        list_options += item + "\n"
    }
    """
    echo "TractoFlow pipeline\n" >> readme.txt
    echo "Start time: $workflow.start\n" >> readme.txt
    echo "[Command-line]\n$workflow.commandLine\n" >> readme.txt
    echo "[Git Info]\n" >> readme.txt
    echo "$workflow.repository - $workflow.revision [$workflow.commitId]\n" >> readme.txt
    echo "[Options]\n" >> readme.txt
    echo "$list_options" >> readme.txt
    """
}

process Bet_Prelim_DWI {
    cpus 2

    input:
    set sid, file(dwi), val(rev), file(bval), file(bvec) from dwi_gradient_for_prelim_bet
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__b0_bet_mask_dilated.nii.gz" into\
        b0_mask_for_eddy
    file "${sid}__b0_bet.nii.gz"
    file "${sid}__b0_bet_mask.nii.gz"

    when:
    rev_b0_count == 0 || (!params.run_topup && params.run_eddy)

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0 --force_b0_threshold
    bet ${sid}__b0.nii.gz ${sid}__b0_bet.nii.gz -m -R -f $params.bet_prelim_f
    scil_image_math.py convert ${sid}__b0_bet_mask.nii.gz ${sid}__b0_bet_mask.nii.gz --data_type uint8 -f
    maskfilter ${sid}__b0_bet_mask.nii.gz dilate ${sid}__b0_bet_mask_dilated.nii.gz\
        --npass $params.dilate_b0_mask_prelim_brain_extraction -nthreads 1
    mrcalc ${sid}__b0.nii.gz ${sid}__b0_bet_mask_dilated.nii.gz\
        -mult ${sid}__b0_bet.nii.gz -quiet -force -nthreads 1
    """
}

process Denoise_DWI {
    cpus params.processes_denoise_dwi
    label 'big_mem'

    input:
    set sid, file(dwi), val(rev) from dwi_for_denoise

    output:
    set sid, "${sid}_${rev}dwi_denoised.nii.gz", val(rev) into\
        dwi_denoised_for_mix

    when:
    params.run_dwi_denoising

    script:
    // The denoised DWI is clipped to 0 since negative values
    // could have been introduced.
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    dwidenoise $dwi dwi_denoised.nii.gz -extent $params.extent -nthreads $task.cpus
    fslmaths dwi_denoised.nii.gz -thr 0 ${sid}_${rev}dwi_denoised.nii.gz
    """
}

dwi_for_test_denoise
    .map{it -> if(!params.run_dwi_denoising){it}}
    .mix(dwi_denoised_for_mix)
    .into{dwi_for_gibbs; dwi_for_test_gibbs}

process Gibbs_correction {
    cpus params.processes_denoise_dwi

    input:
    set sid, file(dwi), val(rev) from dwi_for_gibbs

    output:
    set sid, "${sid}_${rev}dwi_gibbs_corrected.nii.gz", val(rev) into\
        dwi_gibbs_for_mix

    when:
        params.run_gibbs_correction

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    mrdegibbs $dwi ${sid}_${rev}dwi_gibbs_corrected.nii.gz -nthreads $task.cpus
    """
}

dwi_for_test_gibbs
    .map{it -> if(!params.run_gibbs_correction){it}}
    .mix(dwi_gibbs_for_mix)
    .into{dwi_for_eddy; dwi_for_topup; dwi_for_eddy_topup; dwi_for_test_eddy_topup}

dwi_for_topup
    .join(gradients_for_prepare_topup)
    .join(rev_b0_for_prepare_topup)
    .set{dwi_gradients_rev_b0_for_prepare_topup}

process Prepare_for_Topup {
  cpus 2

  input:
    set sid, file(dwi), val(rev), file(bval), file(bvec), file(rev_b0)\
      from dwi_gradients_rev_b0_for_prepare_topup

  output:
    set sid, "${rev_b0}", "${sid}__b0_mean.nii.gz" into simple_rev_b0_for_topup

  when:
    params.run_topup && params.run_eddy

  script:
  """
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_mean.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0 --force_b0_threshold
  """
}

simple_rev_b0_for_topup.mix(complex_rev_b0_for_topup).set{rev_b0_for_topup}

rev_b0_for_topup
  .join(readout_encoding_for_topup)
  .set{rev_b0_with_readout_encoding_for_topup}

process Topup {
    cpus 4

    input:
      set sid, file(rev_b0), file(b0),  readout, encoding\
        from rev_b0_with_readout_encoding_for_topup

    output:
      set sid, "${sid}__corrected_b0s.nii.gz", "${params.prefix_topup}_fieldcoef.nii.gz",
      "${params.prefix_topup}_movpar.txt" into topup_files_for_eddy_topup
      file "${sid}__rev_b0_warped.nii.gz"

    when:
      params.run_topup && params.run_eddy

    script:
    """
      export OMP_NUM_THREADS=$task.cpus
      export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
      export OPENBLAS_NUM_THREADS=1
      export ANTS_RANDOM_SEED=1234
      antsRegistrationSyNQuick.sh -d 3 -f $b0 -m $rev_b0 -o output -t r -e 1
      mv outputWarped.nii.gz ${sid}__rev_b0_warped.nii.gz
      scil_prepare_topup_command.py $b0 ${sid}__rev_b0_warped.nii.gz\
          --config $params.config_topup\
          --encoding_direction $encoding\
          --readout $readout --out_prefix $params.prefix_topup\
          --out_script
      sh topup.sh
      cp corrected_b0s.nii.gz ${sid}__corrected_b0s.nii.gz
    """
}

dwi_for_eddy_topup.into{complex_dwi_for_eddy_topup;simple_dwi_for_eddy_topup}

// DATA FOR CONCATENATE
complex_dwi_for_eddy_topup.combine(sid_rev_dwi_included.collect().toList())
        .filter{ (it[0] in it[3]) }
        .map{it -> [it[0], it[1]]}
        .set{dwi_for_prepare_for_eddy}

// DATA DO NOT CONCATENATE
expl1 = Channel.value(0)

simple_dwi_for_eddy_topup.combine(sid_rev_dwi_excluded.collect().toList())
        .filter{ !(it[0] in it[2]) }
        .map{it -> [it[0], it[1]]}
        .join(gradients_for_eddy_topup)
	.merge(expl1)
        .set{simple_dwi_gradients_for_eddy}


dwi_for_prepare_for_eddy.join(gradients_for_prepare_dwi_for_eddy)
  .groupTuple(by:[0])
  .map{it.flatten().toList()}
  .set{dwi_rev_dwi_gradients_for_prepare_dwi_for_eddy}


process Prepare_dwi_for_eddy {
  cpus 2

  input:
    set sid, file(dwi), file(rev_dwi), file(bval), file(rev_bval), file(bvec), \
        file(rev_bvec) from dwi_rev_dwi_gradients_for_prepare_dwi_for_eddy

  output:
    set sid, "${sid}__concatenated_dwi.nii.gz", "${sid}__concatenated_dwi.bval", "${sid}__concatenated_dwi.bvec", env(rev_number_dir) into concatenated_dwi_for_eddy

  when:
    params.run_topup && params.run_eddy

  script:
  """
    scil_concatenate_dwi.py ${sid}__concatenated_dwi.nii.gz ${sid}__concatenated_dwi.bval ${sid}__concatenated_dwi.bvec -f\
      --in_dwis ${dwi} ${rev_dwi} --in_bvals ${bval} ${rev_bval}\
      --in_bvecs ${bvec} ${rev_bvec}

    rev_number_dir=\$(scil_print_header.py ${rev_dwi} --key dim | sed "s/  / /g" | sed "s/  / /g" | rev | cut -d' ' -f4-4 | rev)
  """
}

simple_dwi_gradients_for_eddy.mix(concatenated_dwi_for_eddy)
    .join(topup_files_for_eddy_topup)
    .join(readout_encoding_for_eddy_topup)
    .set{dwi_gradients_mask_topup_files_for_eddy_topup}

process Eddy_Topup {
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), val(number_rev_dwi), file(b0s_corrected),
        file(field), file(movpar), readout, encoding\
        from dwi_gradients_mask_topup_files_for_eddy_topup
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_from_eddy_topup
    set sid, "${sid}__bval_eddy", "${sid}__dwi_eddy_corrected.bvec" into\
        gradients_from_eddy_topup
    file "${sid}__b0_bet_mask.nii.gz"

    when:
    rev_b0_count > 0 && params.run_topup && params.run_eddy

    // Corrected DWI is clipped to ensure there are no negative values
    // introduced by Eddy.
    script:
        slice_drop_flag=""
        if (params.use_slice_drop_correction)
            slice_drop_flag="--slice_drop_correction"
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        mrconvert $b0s_corrected b0_corrected.nii.gz -coord 3 0 -axes 0,1,2 -nthreads 1
        bet b0_corrected.nii.gz ${sid}__b0_bet.nii.gz -m -R\
            -f $params.bet_topup_before_eddy_f
        scil_prepare_eddy_command.py $dwi $bval $bvec ${sid}__b0_bet_mask.nii.gz\
            --topup $params.prefix_topup --eddy_cmd $params.eddy_cmd\
            --b0_thr $params.b0_thr_extract_b0\
            --encoding_direction $encoding\
            --readout $readout --out_script --fix_seed\
	    --n_reverse ${number_rev_dwi}\
            $slice_drop_flag
        sh eddy.sh
        fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
        mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
}

dwi_for_eddy
    .join(gradients_for_eddy)
    .join(b0_mask_for_eddy)
    .join(readout_encoding_for_eddy)
    .set{dwi_gradients_mask_topup_files_for_eddy}

process Eddy {
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), file(mask), readout, encoding\
        from dwi_gradients_mask_topup_files_for_eddy
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_from_eddy
    set sid, "${sid}__bval_eddy", "${sid}__dwi_eddy_corrected.bvec" into\
        gradients_from_eddy

    when:
    (rev_b0_count == 0 && params.run_eddy) || (!params.run_topup && params.run_eddy)

    // Corrected DWI is clipped to 0 since Eddy can introduce negative values.
    script:
        slice_drop_flag=""
        if (params.use_slice_drop_correction) {
            slice_drop_flag="--slice_drop_correction"
        }
        """
        export OMP_NUM_THREADS=$task.cpus
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        export OPENBLAS_NUM_THREADS=1
        scil_prepare_eddy_command.py $dwi $bval $bvec $mask\
            --eddy_cmd $params.eddy_cmd --b0_thr $params.b0_thr_extract_b0\
            --encoding_direction $encoding\
            --readout $readout --out_script --fix_seed\
            $slice_drop_flag
        sh eddy.sh
        fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
        mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
}


dwi_for_test_eddy_topup.map{it -> if(!params.run_eddy){it}}.set{dwi_for_skip_eddy_topup}
gradients_for_test_eddy_topup.map{it -> if(!params.run_eddy){it}}.set{gradients_for_skip_eddy_topup}

dwi_from_eddy
    .mix(dwi_from_eddy_topup)
    .mix(dwi_for_skip_eddy_topup)
    .set{dwi_for_bet}

gradients_from_eddy
    .mix(gradients_from_eddy_topup)
    .mix(gradients_for_skip_eddy_topup)
    .into{gradients_for_extract_b0;
          gradients_for_dti_shell;
          gradients_for_fodf_shell;
          gradients_for_normalize;
          gradients_for_bet;
          gradients_for_sh_fitting_shell}

dwi_for_bet
    .join(gradients_for_bet)
    .set{dwi_gradients_for_bet}

process Bet_DWI {
    cpus 2
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradients_for_bet

    output:
    set sid, "${sid}__b0_bet.nii.gz", "${sid}__b0_bet_mask.nii.gz" into\
        b0_and_mask_for_crop
    set sid, "${sid}__dwi_bet.nii.gz", "${sid}__b0_bet.nii.gz",
        "${sid}__b0_bet_mask.nii.gz" into dwi_b0_b0_mask_for_n4
    file "${sid}__b0_no_bet.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_no_bet.nii.gz --mean\
            --b0_thr $params.b0_thr_extract_b0 --force_b0_threshold
    bet ${sid}__b0_no_bet.nii.gz ${sid}__b0_bet.nii.gz -m -R -f $params.bet_dwi_final_f
    scil_image_math.py convert ${sid}__b0_bet_mask.nii.gz ${sid}__b0_bet_mask.nii.gz --data_type uint8 -f
    mrcalc $dwi ${sid}__b0_bet_mask.nii.gz -mult ${sid}__dwi_bet.nii.gz -quiet -nthreads 1
    """
}

process N4_DWI {
    cpus 1
    label 'big_mem'

    input:
    set sid, file(dwi), file(b0), file(b0_mask)\
        from dwi_b0_b0_mask_for_n4

    output:
    set sid, "${sid}__dwi_n4.nii.gz" into dwi_for_crop

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    N4BiasFieldCorrection -i $b0\
        -o [${sid}__b0_n4.nii.gz, bias_field_b0.nii.gz]\
        -c [300x150x75x50, 1e-6] -v 1
    scil_apply_bias_field_on_dwi.py $dwi bias_field_b0.nii.gz\
        ${sid}__dwi_n4.nii.gz --mask $b0_mask -f
    """
}

dwi_for_crop
    .join(b0_and_mask_for_crop)
    .set{dwi_and_b0_mask_b0_for_crop}

process Crop_DWI {
    cpus 1
    label 'big_mem'

    input:
    set sid, file(dwi), file(b0), file(b0_mask) from dwi_and_b0_mask_b0_for_crop

    output:
    set sid, "${sid}__dwi_cropped.nii.gz",
        "${sid}__b0_mask_cropped.nii.gz" into dwi_mask_for_normalize
    set sid, "${sid}__b0_mask_cropped.nii.gz" into mask_for_resample
    file "${sid}__b0_cropped.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_crop_volume.py $dwi ${sid}__dwi_cropped.nii.gz -f\
        --output_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0 ${sid}__b0_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0_mask ${sid}__b0_mask_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    scil_image_math.py convert ${sid}__b0_mask_cropped.nii.gz ${sid}__b0_mask_cropped.nii.gz --data_type uint8 -f
    """
}

process Denoise_T1 {
    cpus params.processes_denoise_t1

    input:
    set sid, file(t1) from t1_for_denoise

    output:
    set sid, "${sid}__t1_denoised.nii.gz" into t1_for_mix_n4

    when:
    params.run_t1_denoising

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_run_nlmeans.py $t1 ${sid}__t1_denoised.nii.gz 1 \
        --processes $task.cpus -f
    """
}

t1_for_test_denoise
    .map{it -> if(!params.run_t1_denoising){it}}
    .mix(t1_for_mix_n4)
    .set{t1_for_n4}

process N4_T1 {
    cpus 1

    input:
    set sid, file(t1) from t1_for_n4

    output:
    set sid, "${sid}__t1_n4.nii.gz" into t1_for_resample, t1_for_test_resample

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    N4BiasFieldCorrection -i $t1\
        -o [${sid}__t1_n4.nii.gz, bias_field_t1.nii.gz]\
        -c [300x150x75x50, 1e-6] -v 1
    """
}

process Resample_T1 {
    cpus 1

    input:
    set sid, file(t1) from t1_for_resample

    output:
    set sid, "${sid}__t1_resampled.nii.gz" into t1_resampled_for_mix

    when:
    params.run_resample_t1

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_resample_volume.py $t1 ${sid}__t1_resampled.nii.gz \
        --voxel_size $params.t1_resolution \
        --interp  $params.t1_interpolation
    """
}

t1_for_test_resample
    .map{it -> if(!params.run_resample_t1){it}}
    .mix(t1_resampled_for_mix)
    .set{t1_for_bet}

process Bet_T1 {
    cpus params.processes_brain_extraction_t1

    input:
    set sid, file(t1) from t1_for_bet

    output:
    set sid, "${sid}__t1_bet.nii.gz", "${sid}__t1_bet_mask.nii.gz"\
        into t1_and_mask_for_crop

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsBrainExtraction.sh -d 3 -a $t1 -e $params.template_t1/t1_template.nii.gz\
        -o bet/ -m $params.template_t1/t1_brain_probability_map.nii.gz -u 0
    scil_image_math.py convert bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz --data_type uint8
    mrcalc $t1 ${sid}__t1_bet_mask.nii.gz -mult ${sid}__t1_bet.nii.gz -nthreads 1
    """
}

process Crop_T1 {
    cpus 1

    input:
    set sid, file(t1), file(t1_mask) from t1_and_mask_for_crop

    output:
    set sid, "${sid}__t1_bet_cropped.nii.gz", "${sid}__t1_bet_mask_cropped.nii.gz"\
        into t1_and_mask_for_reg

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_crop_volume.py $t1 ${sid}__t1_bet_cropped.nii.gz\
        --output_bbox t1_boundingBox.pkl -f
    scil_crop_volume.py $t1_mask ${sid}__t1_bet_mask_cropped.nii.gz\
        --input_bbox t1_boundingBox.pkl -f
    scil_image_math.py convert ${sid}__t1_bet_mask_cropped.nii.gz ${sid}__t1_bet_mask_cropped.nii.gz --data_type uint8 -f
    """
}


dwi_mask_for_normalize
    .join(gradients_for_normalize)
    .set{dwi_mask_grad_for_normalize}
process Normalize_DWI {
    cpus 3
    label 'big_mem'

    input:
    set sid, file(dwi), file(mask), file(bval), file(bvec) from dwi_mask_grad_for_normalize

    output:
    set sid, "${sid}__dwi_normalized.nii.gz" into dwi_for_resample, dwi_for_test_resample
    file "${sid}_fa_wm_mask.nii.gz"

    script:
    if (params.dti_shells)
      """
      export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
      export OMP_NUM_THREADS=1
      export OPENBLAS_NUM_THREADS=1
      scil_extract_dwi_shell.py $dwi \
          $bval $bvec $params.dti_shells dwi_dti.nii.gz \
          bval_dti bvec_dti -t $params.dwi_shell_tolerance
      scil_compute_dti_metrics.py dwi_dti.nii.gz bval_dti bvec_dti --mask $mask\
          --not_all --fa fa.nii.gz --force_b0_threshold
      mrthreshold fa.nii.gz ${sid}_fa_wm_mask.nii.gz -abs $params.fa_mask_threshold -nthreads 1
      dwinormalise $dwi ${sid}_fa_wm_mask.nii.gz ${sid}__dwi_normalized.nii.gz\
          -fslgrad $bvec $bval -nthreads 1
      """
    else
      """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        shells=\$(cut -d ' ' --output-delimiter=\$'\\n' -f 1- $bval | awk -F' ' '{v=int(\$1)}{if(v<=$params.max_dti_shell_value)print v}' | uniq)

        scil_extract_dwi_shell.py $dwi \
            $bval $bvec \$shells dwi_dti.nii.gz \
            bval_dti bvec_dti -t $params.dwi_shell_tolerance
        scil_compute_dti_metrics.py dwi_dti.nii.gz bval_dti bvec_dti --mask $mask\
            --not_all --fa fa.nii.gz --force_b0_threshold
        mrthreshold fa.nii.gz ${sid}_fa_wm_mask.nii.gz -abs $params.fa_mask_threshold -nthreads 1
        dwinormalise $dwi ${sid}_fa_wm_mask.nii.gz ${sid}__dwi_normalized.nii.gz\
            -fslgrad $bvec $bval -nthreads 1
      """

}

dwi_for_resample
    .join(mask_for_resample)
    .set{dwi_mask_for_resample}
process Resample_DWI {
    cpus 3

    input:
    set sid, file(dwi), file(mask) from dwi_mask_for_resample

    output:
    set sid, "${sid}__dwi_resampled.nii.gz" into\
        dwi_resampled_for_mix

    when:
    params.run_resample_dwi

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_resample_volume.py $dwi \
        dwi_resample.nii.gz \
        --voxel_size $params.dwi_resolution \
        --interp  $params.dwi_interpolation
    fslmaths dwi_resample.nii.gz -thr 0 dwi_resample_clipped.nii.gz
    scil_resample_volume.py $mask \
        mask_resample.nii.gz \
        --ref dwi_resample.nii.gz \
        --enforce_dimensions \
        --interp nn
    mrcalc dwi_resample_clipped.nii.gz mask_resample.nii.gz\
        -mult ${sid}__dwi_resampled.nii.gz -quiet -nthreads 1
    """
}

dwi_for_test_resample
    .map{it -> if(!params.run_resample_dwi){it}}
    .mix(dwi_resampled_for_mix)
    .into{dwi_for_extract_b0; dwi_for_extract_dti_shell; dwi_for_extract_fodf_shell; dwi_for_extract_sh_fitting_shell}

dwi_for_extract_b0
    .join(gradients_for_extract_b0)
    .set{dwi_and_grad_for_extract_b0}

process Extract_B0 {
    cpus 3
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_and_grad_for_extract_b0

    output:
    set sid, "${sid}__b0_resampled.nii.gz" into b0_for_reg
    set sid, "${sid}__b0_mask_resampled.nii.gz" into\
        b0_mask_for_dti_metrics,
        b0_mask_for_fodf,
        b0_mask_for_rf

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_resampled.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0 --force_b0_threshold
    mrthreshold ${sid}__b0_resampled.nii.gz ${sid}__b0_mask_resampled.nii.gz\
        --abs 0.00001 -nthreads 1
    """
}

dwi_for_extract_sh_fitting_shell
    .join(gradients_for_sh_fitting_shell)
    .set{dwi_and_grad_for_extract_sh_fitting_shell}

process Extract_SH_Fitting_Shell {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_sh_fitting_shell

    output:
    set sid, "${sid}__dwi_sh_fitting.nii.gz", "${sid}__bval_sh_fitting",
        "${sid}__bvec_sh_fitting" into \
        dwi_and_grad_for_sh_fitting

    when:
    params.sh_fitting

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.sh_fitting_shells ${sid}__dwi_sh_fitting.nii.gz \
        ${sid}__bval_sh_fitting ${sid}__bvec_sh_fitting -t $params.dwi_shell_tolerance -f
    """
}

process SH_Fitting {
    cpus 1

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_and_grad_for_sh_fitting

    output:
    file "${sid}__dwi_sh.nii.gz"

    when:
    params.sh_fitting

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_sh_from_signal.py --sh_order $params.sh_fitting_order --sh_basis $params.sh_fitting_basis $dwi $bval $bvec ${sid}__dwi_sh.nii.gz
    """
}

dwi_for_extract_dti_shell
    .join(gradients_for_dti_shell)
    .set{dwi_and_grad_for_extract_dti_shell}

process Extract_DTI_Shell {
    cpus 3
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_dti_shell

    output:
    set sid, "${sid}__dwi_dti.nii.gz", "${sid}__bval_dti",
        "${sid}__bvec_dti" into \
        dwi_and_grad_for_dti_metrics, \
        dwi_and_grad_for_rf

    script:
    if (params.dti_shells)
      """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_extract_dwi_shell.py $dwi \
          $bval $bvec $params.dti_shells ${sid}__dwi_dti.nii.gz \
          ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
      """
    else
      """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1

        shells=\$(cut -d ' ' --output-delimiter=\$'\\n' -f 1- $bval | \
                awk -F' ' '{v=int(\$1)}{if(v<=$params.max_dti_shell_value)print v}' | uniq)

        scil_extract_dwi_shell.py $dwi \
          $bval $bvec \$shells ${sid}__dwi_dti.nii.gz \
          ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
      """
}

dwi_and_grad_for_dti_metrics
    .join(b0_mask_for_dti_metrics)
    .set{dwi_and_grad_for_dti_metrics}

process DTI_Metrics {
    cpus 3
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask)\
        from dwi_and_grad_for_dti_metrics

    output:
    file "${sid}__ad.nii.gz"
    file "${sid}__evecs.nii.gz"
    file "${sid}__evecs_v1.nii.gz"
    file "${sid}__evecs_v2.nii.gz"
    file "${sid}__evecs_v3.nii.gz"
    file "${sid}__evals.nii.gz"
    file "${sid}__evals_e1.nii.gz"
    file "${sid}__evals_e2.nii.gz"
    file "${sid}__evals_e3.nii.gz"
    file "${sid}__fa.nii.gz"
    file "${sid}__ga.nii.gz"
    file "${sid}__rgb.nii.gz"
    file "${sid}__md.nii.gz"
    file "${sid}__mode.nii.gz"
    file "${sid}__norm.nii.gz"
    file "${sid}__rd.nii.gz"
    file "${sid}__tensor.nii.gz"
    file "${sid}__nonphysical.nii.gz"
    file "${sid}__pulsation_std_dwi.nii.gz"
    file "${sid}__residual.nii.gz"
    file "${sid}__residual_iqr_residuals.npy"
    file "${sid}__residual_mean_residuals.npy"
    file "${sid}__residual_q1_residuals.npy"
    file "${sid}__residual_q3_residuals.npy"
    file "${sid}__residual_residuals_stats.png"
    file "${sid}__residual_std_residuals.npy"
    set sid, "${sid}__fa.nii.gz", "${sid}__md.nii.gz" into fa_md_for_fodf
    set sid, "${sid}__fa.nii.gz" into\
        fa_for_reg, fa_for_pft_tracking, fa_for_local_tracking_mask, fa_for_local_seeding_mask

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_dti_metrics.py $dwi $bval $bvec --mask $b0_mask\
        --ad ${sid}__ad.nii.gz --evecs ${sid}__evecs.nii.gz\
        --evals ${sid}__evals.nii.gz --fa ${sid}__fa.nii.gz\
        --ga ${sid}__ga.nii.gz --rgb ${sid}__rgb.nii.gz\
        --md ${sid}__md.nii.gz --mode ${sid}__mode.nii.gz\
        --norm ${sid}__norm.nii.gz --rd ${sid}__rd.nii.gz\
        --tensor ${sid}__tensor.nii.gz\
        --non-physical ${sid}__nonphysical.nii.gz\
        --pulsation ${sid}__pulsation.nii.gz\
        --residual ${sid}__residual.nii.gz\
        -f --force_b0_threshold
    """
}

dwi_for_extract_fodf_shell
    .join(gradients_for_fodf_shell)
    .set{dwi_and_grad_for_extract_fodf_shell}

process Extract_FODF_Shell {
    cpus 3
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_fodf_shell

    output:
    set sid, "${sid}__dwi_fodf.nii.gz", "${sid}__bval_fodf",
        "${sid}__bvec_fodf" into\
        dwi_and_grad_for_fodf

    script:
    if (params.fodf_shells)
      """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_extract_dwi_shell.py $dwi \
          $bval $bvec $params.fodf_shells ${sid}__dwi_fodf.nii.gz \
          ${sid}__bval_fodf ${sid}__bvec_fodf -t $params.dwi_shell_tolerance -f
      """
    else
      """
      export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
      export OMP_NUM_THREADS=1
      export OPENBLAS_NUM_THREADS=1

      shells=\$(cut -d ' ' --output-delimiter=\$'\\n' -f 1- $bval | \
      awk -F' ' '{v=int(\$1)}{if(v>=$params.min_fodf_shell_value|| \
      v<=$params.b0_thr_extract_b0)print v}' | uniq)

      scil_extract_dwi_shell.py $dwi \
        $bval $bvec \$shells ${sid}__dwi_fodf.nii.gz \
        ${sid}__bval_fodf ${sid}__bvec_fodf -t $params.dwi_shell_tolerance -f
      """
}

t1_and_mask_for_reg
    .join(fa_for_reg)
    .join(b0_for_reg)
    .set{t1_fa_b0_for_reg}

process Register_T1 {
    cpus params.processes_registration

    input:
    set sid, file(t1), file(t1_mask), file(fa), file(b0) from t1_fa_b0_for_reg

    output:
    set sid, "${sid}__t1_warped.nii.gz" into t1_for_seg
    set sid, "${sid}__t1_warped.nii.gz", "${sid}__output0GenericAffine.mat",
        "${sid}__output1Warp.nii.gz" into t1_for_freesurfer_reg
    file "${sid}__output1InverseWarp.nii.gz"
    file "${sid}__t1_mask_warped.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsRegistration --dimensionality 3 --float 0\
        --output [output,outputWarped.nii.gz,outputInverseWarped.nii.gz]\
        --interpolation Linear --use-histogram-matching 0\
        --winsorize-image-intensities [0.005,0.995]\
        --initial-moving-transform [$b0,$t1,1]\
        --transform Rigid['0.2']\
        --metric MI[$b0,$t1,1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform Affine['0.2']\
        --metric MI[$b0,$t1,1,32,Regular,0.25]\
        --convergence [500x250x125x50,1e-6,10] --shrink-factors 8x4x2x1\
        --smoothing-sigmas 3x2x1x0\
        --transform SyN[0.1,3,0]\
        --metric MI[$b0,$t1,1,32]\
        --metric CC[$fa,$t1,1,4]\
        --convergence [50x25x10,1e-6,10] --shrink-factors 4x2x1\
        --smoothing-sigmas 3x2x1
    mv outputWarped.nii.gz ${sid}__t1_warped.nii.gz
    mv output0GenericAffine.mat ${sid}__output0GenericAffine.mat
    mv output1InverseWarp.nii.gz ${sid}__output1InverseWarp.nii.gz
    mv output1Warp.nii.gz ${sid}__output1Warp.nii.gz
    antsApplyTransforms -d 3 -i $t1_mask -r ${sid}__t1_warped.nii.gz \
        -o ${sid}__t1_mask_warped.nii.gz -n NearestNeighbor \
        -t ${sid}__output1Warp.nii.gz ${sid}__output0GenericAffine.mat
    scil_image_math.py convert ${sid}__t1_mask_warped.nii.gz ${sid}__t1_mask_warped.nii.gz --data_type uint8 -f
    """
}

labels_for_reg
    .join(t1_for_freesurfer_reg)
    .set{labels_mat_for_reg}

process Register_Freesurfer {
    cpus 1

    input:
    set sid, file(aparc), file(wmparc), file(t1), file(affine),
        file(warp) from labels_mat_for_reg

    output:
    set sid, "${sid}__aparc_warped.nii.gz", "${sid}__wmparc_warped.nii.gz" \
        into labels_for_segmentation

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export ANTS_RANDOM_SEED=1234
    antsApplyTransforms -d 3 -i $aparc -r ${sid}__t1_warped.nii.gz \
        -o ${sid}__aparc_warped.nii.gz -n NearestNeighbor \
        -t ${warp} ${affine}
    antsApplyTransforms -d 3 -i $wmparc -r ${sid}__t1_warped.nii.gz \
        -o ${sid}__wmparc_warped.nii.gz -n NearestNeighbor \
        -t ${warp} ${affine}
    """
}

process Segment_Freesurfer {
    cpus 1

    input:
    set sid, file(aparc), file(wmparc) from labels_for_segmentation

    output:
    set sid, "${sid}__mask_wm.nii.gz" into wm_mask_freesurfer
    file "${sid}__mask_gm.nii.gz"
    file "${sid}__mask_csf.nii.gz"

    when:
        params.run_tractoflow_abs

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    mkdir wmparc_desikan/
    mkdir wmparc_subcortical/
    mkdir aparc+aseg_subcortical/
    scil_image_math.py convert $aparc aparc+aseg_int16.nii.gz --data_type int16 -f
    scil_image_math.py convert $wmparc wmparc_int16.nii.gz --data_type int16 -f

    scil_split_volume_by_labels.py wmparc_int16.nii.gz --scilpy_lut freesurfer_desikan_killiany --out_dir wmparc_desikan
    scil_split_volume_by_labels.py wmparc_int16.nii.gz --scilpy_lut freesurfer_subcortical --out_dir wmparc_subcortical
    scil_split_volume_by_labels.py aparc+aseg_int16.nii.gz --scilpy_lut freesurfer_subcortical --out_dir aparc+aseg_subcortical

    scil_image_math.py union wmparc_desikan/*\
                             wmparc_subcortical/right-cerebellum-cortex.nii.gz\
                             wmparc_subcortical/left-cerebellum-cortex.nii.gz\
                             mask_cortex_m.nii.gz -f
    scil_image_math.py union wmparc_subcortical/corpus-callosum-*\
                             aparc+aseg_subcortical/*white-matter*\
                             wmparc_subcortical/brain-stem.nii.gz\
                             aparc+aseg_subcortical/*ventraldc*\
                             mask_wm_m.nii.gz -f
    scil_image_math.py union wmparc_subcortical/*thalamus*\
                             wmparc_subcortical/*putamen*\
                             wmparc_subcortical/*pallidum*\
                             wmparc_subcortical/*hippocampus*\
                             wmparc_subcortical/*caudate*\
                             wmparc_subcortical/*amygdala*\
                             wmparc_subcortical/*accumbens*\
                             wmparc_subcortical/*plexus*\
                             mask_nuclei_m.nii.gz -f
    scil_image_math.py union wmparc_subcortical/*-lateral-ventricle.nii.gz\
                             wmparc_subcortical/*-inferior-lateral-ventricle.nii.gz\
                             wmparc_subcortical/cerebrospinal-fluid.nii.gz\
                             wmparc_subcortical/*th-ventricle.nii.gz\
                             mask_csf_1_m.nii.gz -f
    scil_image_math.py lower_threshold mask_wm_m.nii.gz 0.1\
                                          ${sid}__mask_wm_bin.nii.gz -f
    scil_image_math.py lower_threshold mask_cortex_m.nii.gz 0.1\
                                          ${sid}__mask_gm.nii.gz -f
    scil_image_math.py lower_threshold mask_nuclei_m.nii.gz 0.1\
                                          ${sid}__mask_nuclei_bin.nii.gz -f
    scil_image_math.py lower_threshold mask_csf_1_m.nii.gz 0.1\
                                          ${sid}__mask_csf.nii.gz -f
    scil_image_math.py addition ${sid}__mask_wm_bin.nii.gz\
                                ${sid}__mask_nuclei_bin.nii.gz\
                                ${sid}__mask_wm.nii.gz --data_type int16

    scil_image_math.py convert ${sid}__mask_wm.nii.gz ${sid}__mask_wm.nii.gz --data_type uint8 -f
    scil_image_math.py convert ${sid}__mask_gm.nii.gz ${sid}__mask_gm.nii.gz --data_type uint8 -f
    scil_image_math.py convert ${sid}__mask_csf.nii.gz ${sid}__mask_csf.nii.gz --data_type uint8 -f
    """
}

process Segment_Tissues {
    cpus 1

    input:
    set sid, file(t1) from t1_for_seg

    output:
    set sid, "${sid}__map_wm.nii.gz", "${sid}__map_gm.nii.gz",
        "${sid}__map_csf.nii.gz" into map_wm_gm_csf_for_pft_maps
    set sid, "${sid}__mask_wm.nii.gz" into wm_mask_for_pft_tracking, wm_mask_fast
    file "${sid}__mask_gm.nii.gz"
    file "${sid}__mask_csf.nii.gz"

    when:
        !params.run_tractoflow_abs

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    fast -t 1 -n $params.number_of_tissues\
         -H 0.1 -I 4 -l 20.0 -g -o t1.nii.gz $t1
    scil_image_math.py convert t1_seg_2.nii.gz ${sid}__mask_wm.nii.gz --data_type uint8
    scil_image_math.py convert t1_seg_1.nii.gz ${sid}__mask_gm.nii.gz --data_type uint8
    scil_image_math.py convert t1_seg_0.nii.gz ${sid}__mask_csf.nii.gz --data_type uint8
    mv t1_pve_2.nii.gz ${sid}__map_wm.nii.gz
    mv t1_pve_1.nii.gz ${sid}__map_gm.nii.gz
    mv t1_pve_0.nii.gz ${sid}__map_csf.nii.gz
    """
}

wm_mask_freesurfer
    .concat(wm_mask_fast)
    .into{wm_mask_for_local_tracking_mask;wm_mask_for_local_seeding_mask}

dwi_and_grad_for_rf
    .join(b0_mask_for_rf)
    .set{dwi_b0_for_rf}

process Compute_FRF {
    cpus 3
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask)\
        from dwi_b0_for_rf

    output:
    set sid, "${sid}__frf.txt" into unique_frf, unique_frf_for_mean
    file "${sid}__frf.txt" into all_frf_to_collect

    script:
    if (params.set_frf)
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_ssst_frf.py $dwi $bval $bvec frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radii $params.roi_radius --force_b0_threshold
        scil_set_response_function.py frf.txt $params.manual_frf ${sid}__frf.txt
        """
    else
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_ssst_frf.py $dwi $bval $bvec ${sid}__frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radii $params.roi_radius --force_b0_threshold
        """
}

all_frf_to_collect
    .collect()
    .set{all_frf_for_mean_frf}

process Mean_FRF {
    cpus 1
    publishDir = params.Mean_FRF_Publish_Dir
    tag = {"All_FRF"}

    input:
    file(all_frf) from all_frf_for_mean_frf

    output:
    file "mean_frf.txt" into mean_frf

    when:
    params.mean_frf && !params.set_frf

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_mean_frf.py $all_frf mean_frf.txt
    """
}

frf_for_fodf = unique_frf

if (params.mean_frf) {
    frf_for_fodf = unique_frf_for_mean
                   .merge(mean_frf)
                   .map{it -> [it[0], it[2]]}
}

dwi_and_grad_for_fodf
    .join(b0_mask_for_fodf)
    .join(fa_md_for_fodf)
    .join(frf_for_fodf)
    .set{dwi_b0_metrics_frf_for_fodf}

process FODF_Metrics {
    cpus params.processes_fodf
    label 'big_mem'

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask), file(fa),
        file(md), file(frf) from dwi_b0_metrics_frf_for_fodf

    output:
    set sid, "${sid}__fodf.nii.gz" into fodf_for_pft_tracking, fodf_for_local_tracking
    file "${sid}__peaks.nii.gz"
    file "${sid}__peak_indices.nii.gz"
    file "${sid}__afd_max.nii.gz"
    file "${sid}__afd_total.nii.gz"
    file "${sid}__afd_sum.nii.gz"
    file "${sid}__nufo.nii.gz"

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_ssst_fodf.py $dwi $bval $bvec $frf ${sid}__fodf.nii.gz\
        --sh_order $params.sh_order --sh_basis $params.basis --force_b0_threshold\
        --mask $b0_mask --processes $task.cpus

    scil_compute_fodf_max_in_ventricles.py ${sid}__fodf.nii.gz $fa $md\
        --max_value_output ventricles_fodf_max_value.txt --sh_basis $params.basis\
        --fa_t $params.max_fa_in_ventricle --md_t $params.min_md_in_ventricle\
        -f

    a_threshold=\$(echo $params.fodf_metrics_a_factor*\$(cat ventricles_fodf_max_value.txt)|bc)

    scil_compute_fodf_metrics.py ${sid}__fodf.nii.gz\
        --mask $b0_mask --sh_basis $params.basis\
        --peaks ${sid}__peaks.nii.gz --peak_indices ${sid}__peak_indices.nii.gz\
        --afd_max ${sid}__afd_max.nii.gz --afd_total ${sid}__afd_total.nii.gz\
        --afd_sum ${sid}__afd_sum.nii.gz --nufo ${sid}__nufo.nii.gz\
        --rt $params.relative_threshold --at \${a_threshold}
    """
}

process PFT_Tracking_Maps {
    cpus 1

    input:
    set sid, file(wm), file(gm), file(csf) from map_wm_gm_csf_for_pft_maps

    output:
    set sid, "${sid}__map_include.nii.gz",
        "${sid}__map_exclude.nii.gz" into pft_maps_for_pft_tracking
    set sid, "${sid}__interface.nii.gz" into interface_for_pft_seeding_mask

    when:
        params.run_pft_tracking

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    scil_compute_maps_for_particle_filter_tracking.py $wm $gm $csf \
        --include ${sid}__map_include.nii.gz \
        --exclude ${sid}__map_exclude.nii.gz \
        --interface ${sid}__interface.nii.gz -f
    """
}
wm_mask_for_pft_tracking
    .join(fa_for_pft_tracking)
    .join(interface_for_pft_seeding_mask)
    .set{wm_fa_int_for_pft}

process PFT_Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa), file(interface_mask) from wm_fa_int_for_pft

    output:
    set sid, "${sid}__pft_seeding_mask.nii.gz" into seeding_mask_for_pft

    when:
        params.run_pft_tracking

    script:
    if (params.pft_seeding_mask_type == "wm")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_image_math.py union $wm $interface_mask ${sid}__pft_seeding_mask.nii.gz\
            --data_type uint8
        """
    else if (params.pft_seeding_mask_type == "interface")
        """
        mv $interface_mask ${sid}__pft_seeding_mask.nii.gz
        """
    else if (params.pft_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.pft_fa_seeding_mask_threshold -ge ${sid}__pft_seeding_mask.nii.gz\
          -datatype uint8
        """
}

fodf_for_pft_tracking
    .join(pft_maps_for_pft_tracking)
    .join(seeding_mask_for_pft)
    .set{fodf_maps_for_pft_tracking}

process PFT_Tracking {
    cpus 2

    input:
    set sid, file(fodf), file(include), file(exclude), file(seed)\
        from fodf_maps_for_pft_tracking
    each curr_seed from pft_random_seed

    output:
    file "${sid}__pft_tracking_${params.pft_algo}_${params.pft_seeding_mask_type}_seed_${curr_seed}.trk"

    when:
        params.run_pft_tracking

    script:
    compress =\
        params.pft_compress_streamlines ? '--compress ' + params.pft_compress_value : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_pft.py $fodf $seed $include $exclude\
            ${sid}__pft_tracking_${params.pft_algo}_${params.pft_seeding_mask_type}_seed_${curr_seed}.trk\
            --algo $params.pft_algo --$params.pft_seeding $params.pft_nbr_seeds\
            --seed $curr_seed --step $params.pft_step --theta $params.pft_theta\
            --sfthres $params.pft_sfthres --sfthres_init $params.pft_sfthres_init\
            --min_length $params.pft_min_len --max_length $params.pft_max_len\
            --particles $params.pft_particles --back $params.pft_back\
            --forward $params.pft_front $compress --sh_basis $params.basis
        """
}

wm_mask_for_local_tracking_mask
    .join(fa_for_local_tracking_mask)
    .set{wm_fa_for_local_tracking_mask}

process Local_Tracking_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa) from wm_fa_for_local_tracking_mask

    output:
    set sid, "${sid}__local_tracking_mask.nii.gz" into tracking_mask_for_local

    when:
        params.run_local_tracking

    script:
    if (params.local_tracking_mask_type == "wm")
        """
        mv $wm ${sid}__local_tracking_mask.nii.gz
        """
    else if (params.local_tracking_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_tracking_mask_threshold -ge ${sid}__local_tracking_mask.nii.gz\
          -datatype uint8
        """
}

wm_mask_for_local_seeding_mask
    .join(fa_for_local_seeding_mask)
    .set{wm_fa_for_local_seeding_mask}

process Local_Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(fa) from wm_fa_for_local_seeding_mask

    output:
    set sid, "${sid}__local_seeding_mask.nii.gz" into tracking_seeding_mask_for_local

    when:
        params.run_local_tracking

    script:
    if (params.local_seeding_mask_type == "wm")
        """
        mv $wm ${sid}__local_seeding_mask.nii.gz
        """
    else if (params.local_seeding_mask_type == "fa")
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        mrcalc $fa $params.local_fa_seeding_mask_threshold -ge ${sid}__local_seeding_mask.nii.gz -datatype uint8
        """
}

fodf_for_local_tracking
    .join(tracking_mask_for_local)
    .join(tracking_seeding_mask_for_local)
    .set{fodf_maps_for_local_tracking}

process Local_Tracking {
    cpus 2

    input:
    set sid, file(fodf), file(tracking_mask), file(seed)\
        from fodf_maps_for_local_tracking
    each curr_seed from local_random_seed

    output:
    file "${sid}__local_tracking_${params.local_algo}_${params.local_seeding_mask_type}_seeding_${params.local_tracking_mask_type}_mask_seed_${curr_seed}.trk"

    when:
        params.run_local_tracking

    script:
    compress =\
        params.local_compress_streamlines ? '--compress ' + params.local_compress_value : ''
        """
        export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        scil_compute_local_tracking.py $fodf $seed $tracking_mask\
            ${sid}__local_tracking_${params.local_algo}_${params.local_seeding_mask_type}_seeding_${params.local_tracking_mask_type}_mask_seed_${curr_seed}.trk\
            --algo $params.local_algo --$params.local_seeding $params.local_nbr_seeds\
            --seed $curr_seed --step $params.local_step --theta $params.local_theta\
            --sfthres $params.local_sfthres --min_length $params.local_min_len\
            --max_length $params.local_max_len $compress --sh_basis $params.basis
        """
}
