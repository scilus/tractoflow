#!/usr/bin/env nextflow

import groovy.json.*

params.input = false
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
    .into{rev_b0; check_rev_b0}
}
else if (params.bids || params.bids_config){
    if (!params.bids_config) {
        log.info "Input BIDS: $params.bids"
        if (params.participants_label) {
            Integer start = workflow.commandLine.indexOf("participants_label") + "participants_label".length();
            participant_cleaned = workflow.commandLine.substring(start, workflow.commandLine.indexOf("--", start) == -1 ? workflow.commandLine.length() : workflow.commandLine.indexOf("--", start)).replace("=", "").replace("\'", "")
            log.info "Participants: $participant_cleaned"
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

            clean_flag = params.clean_bids ? '--clean ' : ''

            """
            scil_validate_bids.py $bids_folder tractoflow_bids_struct.json\
                --readout $params.readout $participants_flag $clean_flag
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
            sub = [sid, file(item.bval), file(item.bvec), file(item.dwi),
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

            if(item.rev_dwi){
              ch_rev_in_data = [sid, file(item.rev_bval), file(item.rev_bvec), file(item.rev_dwi),
                                file(item.t1), item.TotalReadoutTime, item.DWIPhaseEncodingDir[0]]
              ch_sid_rev_dwi.bind([sid])
              ch_in_data.bind(ch_rev_in_data)
            }
        }
        ch_sid_rev_dwi.close()
        ch_in_data.close()
        ch_simple_rev_b0.close()
        ch_complex_rev_b0.close()
    }

    ch_sid_rev_dwi.into{sid_rev_dwi_included; sid_rev_dwi_excluded}
    ch_in_data.into{in_data; check_subjects_number; in_data_print}

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

if (params.bids && workflow.profile.contains("ABS")){
    error "Error ~ --bids parameter cannot be run with Atlas Based Segmentation (ABS) profile"
}

(dwi, gradients, t1, readout_encoding) = in_data
    .map{sid, bvals, bvecs, dwi, t1, readout, encoding -> [tuple(sid, dwi),
                                        tuple(sid, bvals, bvecs),
                                        tuple(sid, t1),
                                        tuple(sid, readout, encoding)]}
    .separate(4)

t1.unique()
    .into{t1_for_denoise; t1_for_test_denoise}

check_complex_rev_b0.concat(check_simple_rev_b0).count().into{rev_b0_counter; number_rev_b0_for_compare}

unique_subjects_number.count().into{ number_subj_for_null_check; number_subj_for_compare}

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

dwi.into{dwi_for_prelim_bet; dwi_for_denoise; dwi_for_test_denoise; dwi_for_eddy_topup}

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
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradient_for_prelim_bet
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
    set sid, file(dwi) from dwi_for_denoise

    output:
    set sid, "${sid}__dwi_denoised.nii.gz" into\
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
    fslmaths dwi_denoised.nii.gz -thr 0 ${sid}__dwi_denoised.nii.gz
    """
}

dwi_for_test_denoise
    .map{it -> if(!params.run_dwi_denoising){it}}
    .mix(dwi_denoised_for_mix)
    .into{dwi_for_gibbs; dwi_for_merge_dwi; dwi_for_eddy; dwi_for_topup; dwi_for_eddy_topup; dwi_for_test_gibbs}

process Gibbs_correction {
    cpus params.processes_denoise_dwi

    input:
    set sid, file(dwi) from dwi_for_gibbs

    output:
    set sid, "${sid}__dwi_gibbs_corrected.nii.gz" into\
        dwi_gibbs_for_mix

    when:
        params.run_gibbs_correction

    script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=1
    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    mrdegibbs $dwi ${sid}__dwi_gibbs_corrected.nii.gz -nthreads $task.cpus
    """
}

dwi_for_test_gibbs
    .map{it -> if(!params.run_gibbs_correction){it}}
    .mix(dwi_gibbs_for_mix)
    .into{dwi_for_merge_dwi; dwi_for_eddy; dwi_for_topup; dwi_for_eddy_topup; dwi_for_test_eddy_topup}

dwi_for_topup
    .join(gradients_for_prepare_topup)
    .join(rev_b0_for_prepare_topup)
    .set{dwi_gradients_rev_b0_for_prepare_topup}

process Prepare_for_Topup {
  cpus 2

  input:
    set sid, file(dwi), file(bval), file(bvec), file(rev_b0)\
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
        .filter{ (it[0] in it[2]) }
        .map{it -> [it[0], it[1]]}
        .set{dwi_for_prepare_for_eddy}

// DATA DO NOT CONCATENATE
simple_dwi_for_eddy_topup.combine(sid_rev_dwi_excluded.collect().toList())
        .filter{ !(it[0] in it[2]) }
        .map{it -> [it[0], it[1]]}
        .join(gradients_for_eddy_topup)
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
    set sid, "${sid}__concatenated_dwi.nii.gz", "${sid}__concatenated_dwi.bval", "${sid}__concatenated_dwi.bvec" into concatenated_dwi_for_eddy

  when:
    params.run_topup && params.run_eddy

  script:
  """
    scil_concatenate_dwi.py ${sid}__concatenated_dwi.nii.gz ${sid}__concatenated_dwi.bval ${sid}__concatenated_dwi.bvec -f\
      --in_dwis ${dwi} ${rev_dwi} --in_bvals ${bval} ${rev_bval}\
      --in_bvecs ${bvec} ${rev_bvec}
  """
}

simple_dwi_gradients_for_eddy.mix(concatenated_dwi_for_eddy)
    .join(topup_files_for_eddy_topup)
    .join(readout_encoding_for_eddy_topup)
    .set{dwi_gradients_mask_topup_files_for_eddy_topup}
