#!/usr/bin/env nextflow

params.root = false
params.subject = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["b0_thr_extract_b0":"$params.b0_thr_extract_b0",
                "dwi_shell_tolerance":"$params.dwi_shell_tolerance",
                "dilate_b0_mask_prelim_brain_extraction":"$params.dilate_b0_mask_prelim_brain_extraction",
                "bet_prelim_f":"$params.bet_prelim_f",
                "run_dwi_denoising":"$params.run_dwi_denoising",
                "extent":"$params.extent",
                "run_topup":"$params.run_topup",
                "encoding_direction":"$params.encoding_direction",
                "dwell_time":"$params.dwell_time",
                "run_eddy":"$params.run_eddy",
                "eddy_cmd":"$params.eddy_cmd",
                "bet_topup_before_eddy_f":"$params.bet_topup_before_eddy_f",
                "use_slice_drop_correction":"$params.use_slice_drop_correction",
                "bet_dwi_final_f":"$params.bet_dwi_final_f",
                "run_resample_dwi":"$params.run_resample_dwi",
                "dwi_resolution":"$params.dwi_resolution",
                "dwi_interpolation":"$params.dwi_interpolation",
                "dti_shells":"$params.dti_shells",
                "fodf_shells":"$params.fodf_shells",
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
                "wm_seeding":"$params.wm_seeding",
                "algo":"$params.algo",
                "seeding":"$params.seeding",
                "nbr_seeds":"$params.nbr_seeds",
                "random":"$params.random",
                "step":"$params.step",
                "rk_order":"$params.rk_order",
                "theta":"$params.theta",
                "minL":"$params.minL",
                "maxL":"$params.maxL",
                "compress_streamlines":"$params.compress_streamlines",
                "compress_value":"$params.compress_value"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "SCIL Human processing pipeline"
log.info "=============================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Denoise DWI]"
log.info "Denoise DWI: $params.run_dwi_denoising"
log.info ""
log.info "[Topup]"
log.info "Run Topup: $params.run_topup"
log.info ""
log.info "[Eddy]"
log.info "Run Eddy: $params.run_eddy"
log.info "Eddy command: $params.eddy_cmd"
log.info ""
log.info "[Resample DWI]"
log.info "Resample DWI: $params.run_resample_dwi"
log.info "Resolution: $params.dwi_resolution"
log.info ""
log.info "[DTI shells]"
log.info "DTI shells: $params.dti_shells"
log.info ""
log.info "[fODF shells]"
log.info "fODF shells: $params.fodf_shells"
log.info ""
log.info "[Compute fiber response function (FRF)]"
log.info "Set FRF: $params.set_frf"
log.info "FRF value: $params.manual_frf"
log.info ""
log.info "[Mean FRF]"
log.info "Mean FRF: $params.mean_frf"
log.info ""
log.info "[Seeding mask]"
log.info "WM seeding: $params.wm_seeding"
log.info ""
log.info "[PFT tracking]"
log.info "Algo: $params.algo"
log.info "Seeding type: $params.seeding"
log.info "Number of seeds: $params.nbr_seeds"
log.info "Random seed: $params.random"
log.info "Step size: $params.step"
log.info "Theta: $params.theta"
log.info "Minimum length: $params.minL"
log.info "Maximum length: $params.maxL"
log.info "Compress streamlines: $params.compress_streamlines"
log.info "Compressing threshold: $params.compress_value"
log.info ""

log.info "Number of processes per tasks"
log.info "============================="
log.info "DWI brain extraction: $params.processes_brain_extraction_dwi"
log.info "T1 brain extraction: $params.processes_brain_extraction_t1"
log.info "Denoise DWI: $params.processes_denoise_dwi"
log.info "Denoise T1: $params.processes_denoise_t1"
log.info "Eddy: $params.processes_eddy"
log.info "Compute fODF: $params.processes_fodf"
log.info "Registration: $params.processes_registration"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

Channel
    .fromPath("$params.template_t1", type:'dir')
    .set{template_dir_t1}

if (params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    in_data = Channel
        .fromFilePairs("$root/**/*{bval,bvec,dwi.nii.gz,t1.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}
    Channel
    .fromPath("$root/**/*rev_b0.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{rev_b0; check_rev_b0}
    }
else if (params.subject){
    log.info "Input: $params.subject"
    subject = file(params.subject)
    in_data = Channel
        .fromFilePairs("$subject/*{bval,bvec,dwi.nii.gz,t1.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}
    Channel
    .fromPath("$subject/*rev_b0.nii.gz",
                    maxDepth:1)
    .map{[it.parent.name, it]}
    .into{rev_b0; check_rev_b0}
    }

(dwi, gradients, t1_for_denoise) = in_data
    .map{sid, bvals, bvecs, dwi, t1 -> [tuple(sid, dwi),
                                        tuple(sid, bvals, bvecs),
                                        tuple(sid, t1)]}
    .separate(3)

check_rev_b0.count().set{ rev_b0_counter }

dwi.into{dwi_for_prelim_bet; dwi_for_denoise}

gradients
    .into{gradients_for_prelim_bet; gradients_for_eddy; gradients_for_topup;
          gradients_for_eddy_topup}

dwi_for_prelim_bet
    .join(gradients_for_prelim_bet)
    .set{dwi_gradient_for_prelim_bet}

process Bet_Prelim_DWI {
    cpus params.processes_brain_extraction_dwi

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradient_for_prelim_bet

    output:
    set sid, "${sid}__b0_bet_mask_dilated.nii.gz" into\
        b0_mask_for_eddy
    file "${sid}__b0_bet.nii.gz"
    file "${sid}__b0_bet_mask.nii.gz"

    script:
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    bet ${sid}__b0.nii.gz ${sid}__b0_bet.nii.gz -m -R -f $params.bet_prelim_f
    maskfilter ${sid}__b0_bet_mask.nii.gz dilate ${sid}__b0_bet_mask_dilated.nii.gz\
        --npass $params.dilate_b0_mask_prelim_brain_extraction
    mrcalc ${sid}__b0.nii.gz ${sid}__b0_bet_mask_dilated.nii.gz\
        -mult ${sid}__b0_bet.nii.gz -quiet -force
    """
}

process Denoise_DWI {
    cpus params.processes_denoise_dwi

    input:
    set sid, file(dwi) from dwi_for_denoise

    output:
    set sid, "${sid}__dwi_denoised.nii.gz" into\
        dwi_for_eddy,
        dwi_for_topup,
        dwi_for_eddy_topup

    script:
    // The denoised DWI is clipped to 0 since negative values 
    // could have been introduced.
    if(params.run_dwi_denoising)
        """
        MRTRIX_NTHREADS=$task.cpus
        dwidenoise $dwi dwi_denoised.nii.gz -extent $params.extent
        fslmaths dwi_denoised.nii.gz -thr 0 ${sid}__dwi_denoised.nii.gz
        """
    else
        """
        mv $dwi ${sid}__dwi_denoised.nii.gz
        """
}

dwi_for_topup
    .join(gradients_for_topup)
    .join(rev_b0)
    .set{dwi_gradients_rev_b0_for_topup}

process Topup {
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec), file(rev_b0)\
        from dwi_gradients_rev_b0_for_topup

    output:
    set sid, "${sid}__corrected_b0s.nii.gz", "${params.prefix_topup}_fieldcoef.nii.gz",
    "${params.prefix_topup}_movpar.txt" into topup_files_for_eddy_topup

    when:
    params.run_topup && params.run_eddy

    script:
    """
    OMP_NUM_THREADS=$task.cpus
    scil_prepare_topup_command.py $dwi $bval $bvec $rev_b0\
        --config $params.config_topup --b0_thr $params.b0_thr_extract_b0\
        --encoding_direction $params.encoding_direction\
        --dwell_time $params.dwell_time --output_prefix $params.prefix_topup\
        --output_script
    sh topup.sh
    cp corrected_b0s.nii.gz ${sid}__corrected_b0s.nii.gz
    """
}

dwi_for_eddy
    .join(gradients_for_eddy)
    .join(b0_mask_for_eddy)
    .set{dwi_gradients_mask_topup_files_for_eddy}

process Eddy {
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), file(mask)\
        from dwi_gradients_mask_topup_files_for_eddy
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__dwi_corrected.nii.gz", "${sid}__bval_eddy",
        "${sid}__dwi_eddy_corrected.bvec" into\
        dwi_gradients_from_eddy
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_from_eddy
    set sid, "${sid}__bval_eddy", "${sid}__dwi_eddy_corrected.bvec" into\
        gradients_from_eddy, gradients_for_normalize

    when:
    rev_b0_count == 0 || !params.run_topup || (!params.run_eddy && params.run_topup)

    // Corrected DWI is clipped to 0 since Eddy can introduce negative values.
    script:
    if (params.run_eddy) {
        slice_drop_flag=""
        if (params.use_slice_drop_correction) {
            slice_drop_flag="--slice_drop_correction"
        }
        """
        OMP_NUM_THREADS=$task.cpus
        scil_prepare_eddy_command.py $dwi $bval $bvec $mask\
            --eddy_cmd $params.eddy_cmd --b0_thr $params.b0_thr_extract_b0\
            --encoding_direction $params.encoding_direction\
            --dwell_time $params.dwell_time --output_script --fix_seed\
            $slice_drop_flag
        sh eddy.sh
        fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
        mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
    else {
        """
        mv $dwi ${sid}__dwi_corrected.nii.gz
        mv $bvec ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
}

dwi_for_eddy_topup
    .join(gradients_for_eddy_topup)
    .join(topup_files_for_eddy_topup)
    .set{dwi_gradients_mask_topup_files_for_eddy_topup}

process Eddy_Topup {
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0s_corrected),
        file(field), file(movpar)\
        from dwi_gradients_mask_topup_files_for_eddy_topup
    val(rev_b0_count) from rev_b0_counter

    output:
    set sid, "${sid}__dwi_corrected.nii.gz", "${sid}__bval_eddy",
        "${sid}__dwi_eddy_corrected.bvec" into\
        dwi_gradients_from_eddy_topup
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_from_eddy_topup
    set sid, "${sid}__bval_eddy", "${sid}__dwi_eddy_corrected.bvec" into\
        gradients_from_eddy_topup, gradients_for_normalize
    file "${sid}__b0_bet_mask.nii.gz"

    when:
    rev_b0_count > 0 && params.run_topup

    // Corrected DWI is clipped to ensure there are no negative values
    // introduced by Eddy.
    script:
    if (params.run_eddy) {
        slice_drop_flag=""
        if (params.use_slice_drop_correction)
            slice_drop_flag="--slice_drop_correction"
        """
        OMP_NUM_THREADS=$task.cpus
        ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
        mrconvert $b0s_corrected b0_corrected.nii.gz -coord 3 0 -axes 0,1,2
        bet b0_corrected.nii.gz ${sid}__b0_bet.nii.gz -m -R\
            -f $params.bet_topup_before_eddy_f
        scil_prepare_eddy_command.py $dwi $bval $bvec ${sid}__b0_bet_mask.nii.gz\
            --topup $params.prefix_topup --eddy_cmd $params.eddy_cmd\
            --b0_thr $params.b0_thr_extract_b0\
            --encoding_direction $params.encoding_direction\
            --dwell_time $params.dwell_time --output_script --fix_seed\
            $slice_drop_flag
        sh eddy.sh
        fslmaths dwi_eddy_corrected.nii.gz -thr 0 ${sid}__dwi_corrected.nii.gz
        mv dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
    else {
        """
        mv $dwi ${sid}__dwi_corrected.nii.gz
        mv $bvec ${sid}__dwi_eddy_corrected.bvec
        mv $bval ${sid}__bval_eddy
        """
    }
}

dwi_gradients_from_eddy
    .mix(dwi_gradients_from_eddy_topup)
    .set{dwi_gradients_for_extract_b0}

dwi_from_eddy
    .mix(dwi_from_eddy_topup)
    .set{dwi_for_bet}

gradients_from_eddy
    .mix(gradients_from_eddy_topup)
    .into{gradients_for_resample_b0;
          gradients_for_dti_shell;
          gradients_for_fodf_shell}

process Extract_B0 {
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradients_for_extract_b0

    output:
    set sid, "${sid}__b0.nii.gz" into b0_for_bet

    script:
    """
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    """
}

dwi_for_bet
    .join(b0_for_bet)
    .set{dwi_b0_for_bet}

process Bet_DWI {
    cpus params.processes_brain_extraction_dwi

    input:
    set sid, file(dwi), file(b0) from dwi_b0_for_bet

    output:
    set sid, "${sid}__b0_bet.nii.gz", "${sid}__b0_bet_mask.nii.gz" into\
        b0_and_mask_for_crop
    set sid, "${sid}__dwi_bet.nii.gz", "${sid}__b0_bet.nii.gz", 
        "${sid}__b0_bet_mask.nii.gz" into dwi_b0_b0_mask_for_n4

    script:
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    bet ${sid}__b0.nii.gz ${sid}__b0_bet.nii.gz -m -R -f $params.bet_dwi_final_f
    mrcalc $dwi ${sid}__b0_bet_mask.nii.gz -mult ${sid}__dwi_bet.nii.gz -quiet
    """
}

process N4_DWI {
    cpus 1

    input:
    set sid, file(dwi), file(b0), file(b0_mask)\
        from dwi_b0_b0_mask_for_n4

    output:
    set sid, "${sid}__dwi_n4.nii.gz" into dwi_for_crop

    script:
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
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

    input:
    set sid, file(dwi), file(b0), file(b0_mask) from dwi_and_b0_mask_b0_for_crop

    output:
    set sid, "${sid}__dwi_cropped.nii.gz",
        "${sid}__b0_mask_cropped.nii.gz" into dwi_mask_for_normalize
    set sid, "${sid}__b0_mask_cropped.nii.gz" into mask_for_resample
    file "${sid}__b0_cropped.nii.gz"

    script:
    """
    scil_crop_volume.py $dwi ${sid}__dwi_cropped.nii.gz -f\
        --output_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0 ${sid}__b0_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0_mask ${sid}__b0_mask_cropped.nii.gz\
        --input_bbox dwi_boundingBox.pkl -f
    """
}

process Denoise_T1 {
    cpus params.processes_denoise_t1

    input:
    set sid, file(t1) from t1_for_denoise

    output:
    set sid, "${sid}__t1_denoised.nii.gz" into t1_for_n4

    script:
    """
    scil_run_nlmeans.py $t1 ${sid}__t1_denoised.nii.gz 1 \
        --noise_est basic --processes $task.cpus -f
    """
}

process N4_T1 {
    cpus 1

    input:
    set sid, file(t1) from t1_for_n4

    output:
    set sid, "${sid}__t1_n4.nii.gz" into t1_for_resample

    script:
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
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
    set sid, "${sid}__t1_resampled.nii.gz" into t1_for_bet

    script:
    if(params.run_resample_t1)
        """
        scil_resample_volume.py $t1 ${sid}__t1_resampled.nii.gz \
            --resolution $params.t1_resolution \
            --interp  $params.t1_interpolation
        """
    else
        """
        mv $t1 ${sid}__t1_resampled.nii.gz
        """
}

process Bet_T1 {
    cpus params.processes_brain_extraction_t1

    input:
    set sid, file(t1) from t1_for_bet
    file template_dir from template_dir_t1.first()

    output:
    set sid, "${sid}__t1_bet.nii.gz", "${sid}__t1_bet_mask.nii.gz"\
        into t1_and_mask_for_crop

    script:
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    antsBrainExtraction.sh -d 3 -a $t1 -e $template_dir/t1_template.nii.gz\
        -o bet/ -m $template_dir/t1_brain_probability_map.nii.gz -u 0
    mrcalc $t1 bet/BrainExtractionMask.nii.gz -mult ${sid}__t1_bet.nii.gz
    mv bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz
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
    scil_crop_volume.py $t1 ${sid}__t1_bet_cropped.nii.gz\
        --output_bbox t1_boundingBox.pkl -f
    scil_crop_volume.py $t1 ${sid}__t1_bet_mask_cropped.nii.gz\
        --input_bbox t1_boundingBox.pkl -f
    """
}


dwi_mask_for_normalize
    .join(gradients_for_normalize)
    .set{dwi_mask_grad_for_normalize}
process Normalize_DWI {
    cpus 1

    input:
    set sid, file(dwi), file(mask), file(bval), file(bvec) from dwi_mask_grad_for_normalize

    output:
    set sid, "{sid}__dwi_normalized.nii.gz" into dwi_for_resample

    script:
    """
    dwinormalise $dwi $mask ${sid}__dwi_normalized.nii.gz -fslgrad $bvec $bval
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
        dwi_for_resample_b0,
        dwi_for_extract_dti_shell,
        dwi_for_extract_fodf_shell

    script:
    if (params.run_resample_dwi)
        """
        scil_resample_volume.py $dwi \
            dwi_resample.nii.gz \
            --resolution $params.dwi_resolution \
            --interp  $params.dwi_interpolation
        fslmaths dwi_resample.nii.gz -thr 0 dwi_resample_clipped.nii.gz
        scil_resample_volume.py $mask \
            mask_resample.nii.gz \
            --ref dwi_resample.nii.gz \
            --enforce_dimensions \
            --interp nn
        mrcalc dwi_resample_clipped.nii.gz mask_resample.nii.gz\
            -mult ${sid}__dwi_resampled.nii.gz -quiet
        """
    else
        """
        mv $dwi dwi_resampled.nii.gz
        """
}

dwi_for_resample_b0
    .join(gradients_for_resample_b0)
    .set{dwi_and_grad_for_resample_b0}

process Resample_B0 {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_and_grad_for_resample_b0

    output:
    set sid, "${sid}__b0_resampled.nii.gz" into b0_for_reg
    set sid, "${sid}__b0_mask_resampled.nii.gz" into\
        b0_mask_for_dti_metrics,
        b0_mask_for_fodf,
        b0_mask_for_rf

    script:
    """
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_resampled.nii.gz --mean\
        --b0_thr $params.b0_thr_extract_b0
    mrthreshold ${sid}__b0_resampled.nii.gz ${sid}__b0_mask_resampled.nii.gz\
        --abs 0.00001
    """
}

dwi_for_extract_dti_shell
    .join(gradients_for_dti_shell)
    .set{dwi_and_grad_for_extract_dti_shell}

process Extract_DTI_Shell {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_dti_shell

    output:
    set sid, "${sid}__dwi_dti.nii.gz", "${sid}__bval_dti",
        "${sid}__bvec_dti" into \
        dwi_and_grad_for_dti_metrics, \
        dwi_and_grad_for_rf

    script:
    """
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.dti_shells ${sid}__dwi_dti.nii.gz \
        ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
    """
}

dwi_and_grad_for_dti_metrics
    .join(b0_mask_for_dti_metrics)
    .set{dwi_and_grad_for_dti_metrics}

process DTI_Metrics {
    cpus 3

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
        fa_for_reg

    script:
    """
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
        -f
    """
}

dwi_for_extract_fodf_shell
    .join(gradients_for_fodf_shell)
    .set{dwi_and_grad_for_extract_fodf_shell}

process Extract_FODF_Shell {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec)\
        from dwi_and_grad_for_extract_fodf_shell

    output:
    set sid, "${sid}__dwi_fodf.nii.gz", "${sid}__bval_fodf",
        "${sid}__bvec_fodf" into\
        dwi_and_grad_for_fodf

    script:
    """
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.fodf_shells ${sid}__dwi_fodf.nii.gz \
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
    file "${sid}__output0GenericAffine.mat"
    file "${sid}__output1InverseWarp.nii.gz"
    file "${sid}__output1Warp.nii.gz"
    file "${sid}__t1_mask_warped.nii.gz"

    script:
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
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
    """
}

process Segment_Tissues {
    cpus 1

    input:
    set sid, file(t1) from t1_for_seg

    output:
    set sid, "${sid}__map_wm.nii.gz", "${sid}__map_gm.nii.gz",
        "${sid}__map_csf.nii.gz" into map_wm_gm_csf_for_pft_maps
    set sid, "${sid}__mask_wm.nii.gz" into wm_mask_for_seeding_mask
    file "${sid}__mask_gm.nii.gz"
    file "${sid}__mask_csf.nii.gz"

    script:
    """
    fast -t 1 -n $params.number_of_tissues\
         -H 0.1 -I 4 -l 20.0 -g -o t1.nii.gz $t1
    mv t1_seg_2.nii.gz ${sid}__mask_wm.nii.gz
    mv t1_seg_1.nii.gz ${sid}__mask_gm.nii.gz
    mv t1_seg_0.nii.gz ${sid}__mask_csf.nii.gz
    mv t1_pve_2.nii.gz ${sid}__map_wm.nii.gz
    mv t1_pve_1.nii.gz ${sid}__map_gm.nii.gz
    mv t1_pve_0.nii.gz ${sid}__map_csf.nii.gz
    """
}

dwi_and_grad_for_rf
    .join(b0_mask_for_rf)
    .set{dwi_b0_for_rf}

process Compute_FRF {
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask)\
        from dwi_b0_for_rf

    output:
    set sid, "${sid}__frf.txt" into unique_frf, unique_frf_for_mean
    file "${sid}__frf.txt" into all_frf_to_collect

    script:
    if (params.set_frf)
        """
        scil_compute_ssst_frf.py $dwi $bval $bvec frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radius $params.roi_radius
        scil_set_response_function.py frf.txt $params.manual_frf ${sid}__frf.txt
        """
    else
        """
        scil_compute_ssst_frf.py $dwi $bval $bvec ${sid}__frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radius $params.roi_radius
        """
}

all_frf_to_collect
    .collect()
    .set{all_frf_for_mean_frf}

process Mean_FRF {
    cpus 1
    publishDir = params.MeanFRFPublishDir
    tag = {"All_FRF"}

    input:
    file(all_frf) from all_frf_for_mean_frf

    output:
    file "mean_frf.txt" into mean_frf

    when:
    params.mean_frf

    script:
    """
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

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask), file(fa),
        file(md), file(frf) from dwi_b0_metrics_frf_for_fodf

    output:
    set sid, "${sid}__fodf.nii.gz" into fodf_for_tracking
    file "${sid}__peaks.nii.gz"
    file "${sid}__peak_indices.nii.gz"
    file "${sid}__afd_max.nii.gz"
    file "${sid}__afd_total.nii.gz"
    file "${sid}__afd_sum.nii.gz"
    file "${sid}__nufo.nii.gz"

    script:
    """ 
    scil_compute_fodf.py $dwi $bval $bvec $frf --sh_order $params.sh_order\
        --basis $params.basis --force_b0_threshold --mask $b0_mask\
        --fodf ${sid}__fodf.nii.gz --peaks ${sid}__peaks.nii.gz\
        --peak_indices ${sid}__peak_indices.nii.gz --processes $task.cpus

    scil_compute_fodf_max_in_ventricles.py ${sid}__fodf.nii.gz $fa $md\
        --max_value_output ventricles_fodf_max_value.txt -f

    a_threshold=\$(echo $params.fodf_metrics_a_factor*\$(cat ventricles_fodf_max_value.txt)|bc)

    scil_compute_fodf_metrics.py ${sid}__fodf.nii.gz \${a_threshold}\
        --mask $b0_mask --afd ${sid}__afd_max.nii.gz\
        --afd_total ${sid}__afd_total.nii.gz --afd_sum ${sid}__afd_sum.nii.gz\
        --nufo ${sid}__nufo.nii.gz -f
    """
}

process PFT_Maps {
    cpus 1

    input:
    set sid, file(wm), file(gm), file(csf) from map_wm_gm_csf_for_pft_maps

    output:
    set sid, "${sid}__map_include.nii.gz",
        "${sid}__map_exclude.nii.gz" into pft_maps_for_tracking
    set sid, "${sid}__interface.nii.gz" into interface_for_seeding_mask

    script:
    """
    scil_compute_maps_for_particle_filter_tracking.py $wm $gm $csf \
    --include ${sid}__map_include.nii.gz --exclude ${sid}__map_exclude.nii.gz\
        --interface ${sid}__interface.nii.gz -f
    """
}

wm_mask_for_seeding_mask
    .join(interface_for_seeding_mask)
    .set{wm_interface_for_seeding_mask}
    
process Seeding_Mask {
    cpus 1

    input:
    set sid, file(wm), file(interface_mask) from wm_interface_for_seeding_mask

    output:
    set sid, "${sid}__seeding_mask.nii.gz" into seeding_mask_for_tracking

    script:
    if (params.wm_seeding)
        """
        scil_mask_math.py union $wm $interface_mask ${sid}__seeding_mask.nii.gz
        """
    else
        """
        mv $interface_mask ${sid}__seeding_mask.nii.gz
        """
}

fodf_for_tracking
    .join(pft_maps_for_tracking)
    .join(seeding_mask_for_tracking)
    .set{fodf_maps_for_tracking}

process Tracking {
    cpus 1

    input:
    set sid, file(fodf), file(include), file(exclude), file(seed)\
        from fodf_maps_for_tracking

    output:
    file "${sid}__tracking.trk"

    script:
    compress =\
        params.compress_streamlines ? '--compress ' + params.compress_value : ''
        """
        scil_compute_pft_dipy.py $fodf $seed $include $exclude\
            ${sid}__tracking.trk --algo $params.algo\
            --$params.seeding $params.nbr_seeds --seed $params.random\
            --step $params.step --theta $params.theta\
            --sfthres $params.sfthres --sfthres_init $params.sfthres_init
            --min_length $params.minL --max_length $params.maxL\
            --particles $params.particles --back $params.back\
            --forward $params.front $compress --basis $params.basis
        """
}
