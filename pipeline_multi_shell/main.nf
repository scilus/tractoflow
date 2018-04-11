#!/usr/bin/env nextflow

params.root = false
params.subject = false
params.help = false
template_dir_t1="$workflow.projectDir/../data/mni_152_sym_09c/t1"
template_dir_b0="$workflow.projectDir/../data/mni_152_sym_09c/b0"
config_eddy="$workflow.projectDir/../data/b02b0.cnf"

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Multi-shell pipeline"
log.info "================================"
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.root){
    log.info "Input: $params.root"
    root = file(params.root)
    in_data = Channel
        .fromFilePairs("$root/**/{bval,bvec,dwi.nii.gz,rev_b0.nii.gz,t1.nii.gz}",
                       size: 5,
                       maxDepth:2,
                       flat: true) { it.parent.parent.name + "_-_" + it.parent.name}
                       }
else if (params.subject){
    log.info "Input: $params.subject"
    subject = file(params.subject)
    in_data = Channel
        .fromFilePairs("$subject/{bval,bvec,dwi.nii.gz,rev_b0.nii.gz,t1.nii.gz}",
                       size: 5,
                       maxDepth:2,
                       flat: true) { it.parent.name}
                       }


(dwi, bvals, bvecs, rev_b0, t1_for_resample) = in_data
    .map{sid, bvals, bvecs, dwi, rev_b0, t1 -> [tuple(sid, dwi),
                                                tuple(sid, bvals),
                                                tuple(sid, bvecs),
                                                tuple(sid, rev_b0),
                                                tuple(sid, t1)]}
    .separate(5)


dwi.into{dwi_for_prelim_bet; dwi_for_denoise}
bvecs.into{bvecs_for_prelim_bet; bvecs_for_eddy}
bvals.into{bvals_for_prelim_bet; bvals_for_eddy; bvals_for_extract_b0; bvals_for_resample_b0; bvals_for_extract_dti_shell; bvals_for_extract_fodf_shell}

dwi_for_prelim_bet
    .phase(bvals_for_prelim_bet)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(bvecs_for_prelim_bet)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_gradient_for_prelim_bet}

process bet_prelim_dwi {
    tag { "$sid" }
    cpus params.processes_bet_dwi

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradient_for_prelim_bet

    output:
    set sid, "${sid}__b0_bet_mask_dilate.nii.gz" into b0_mask_for_denoise_dwi
    file "${sid}__b0_bet.nii.gz"
    set sid, "${sid}__b0_bet_mask.nii.gz" into b0_mask_for_eddy

    script:
    dir_id = get_dir(sid)
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --b0_thr $params.b0_thr_extract_b0
    antsBrainExtraction.sh -d 3 -a ${sid}__b0.nii.gz -e $template_dir_b0/b0_template.nii.gz\
        -o bet/ -m $template_dir_b0/b0_brain_probability_map.nii.gz\
        -f $template_dir_b0/b0_brain_registration_mask.nii.gz -k 1
    cp bet/BrainExtractionPriorWarped.nii.gz ${sid}__b0_bet_mask.nii.gz
    maskfilter ${sid}__b0_bet_mask.nii.gz dilate ${sid}__b0_bet_mask_dilate.nii.gz --npass $params.dilate_b0_mask_prelim_bet
    mrcalc ${sid}__b0.nii.gz ${sid}__b0_bet_mask_dilate.nii.gz -mult ${sid}__b0_bet.nii.gz -quiet
    """
}

dwi_for_denoise
    .phase(b0_mask_for_denoise_dwi)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_b0_mask_for_denoise}

process denoise_dwi {
    tag { "$sid" }
    cpus params.processes_denoise_dwi

    input:
    set sid, file(dwi), file(b0_mask) from dwi_b0_mask_for_denoise

    output:
    set sid, "${sid}__dwi_denoised.nii.gz" into dwi_for_eddy

    script:
    dir_id = get_dir(sid)
    if(params.do_denoising_dwi)
        """
        MRTRIX_NTHREADS=$task.cpus
        dwidenoise $dwi ${sid}__dwi_denoised.nii.gz -mask $b0_mask -extent $params.extent
        """
    else
        """
        cp $dwi ${sid}__dwi_denoised.nii.gz
        """
}

dwi_for_eddy
    .phase(bvals_for_eddy)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(bvecs_for_eddy)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(rev_b0)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(b0_mask_for_eddy)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_gradients_rev_b0_mask_for_eddy}

process eddy {
    tag { "$sid" }
    cpus params.processes_eddy

    input:
    set sid, file(dwi), file(bval), file(bvec), file(rev_b0), file(mask) from\
        dwi_gradients_rev_b0_mask_for_eddy

    output:
    set sid, "${sid}__dwi_corrected.nii.gz" into\
        dwi_for_extract_b0, dwi_for_bet
    set sid, "${sid}__dwi_eddy_corrected.bvec" into\
        bvecs_for_extract_b0,
        bvecs_for_denoise,
        bvecs_for_resample_b0,
        bvecs_for_dti_shell,
        bvecs_for_fodf_shell,
        bvecs_for_extract_shell_denoise

    script:
    dir_id = get_dir(sid)
    """
    OMP_NUM_THREADS=$task.cpus
    scil_prepare_eddy_command.py $dwi $bval $bvec $rev_b0 $mask --config $params.config_eddy\
        --eddy_cmd $params.eddy_cmd --b0_thr $params.b0_thr_extract_b0\
        --encoding_direction $params.encoding_direction\
        --dwell_time $params.dwell_time --output_script
    sh eddy.sh
    cp dwi_eddy_corrected.nii.gz ${sid}__dwi_corrected.nii.gz
    cp dwi_eddy_corrected.eddy_rotated_bvecs ${sid}__dwi_eddy_corrected.bvec
    """
}

dwi_for_extract_b0
    .phase(bvals_for_extract_b0)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(bvecs_for_extract_b0)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_gradients_for_extract_b0}

process extract_b0 {
    tag { "$sid" }
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_gradients_for_extract_b0

    output:
    set sid, "${sid}__b0.nii.gz" into b0_for_bet

    script:
    dir_id = get_dir(sid)
    """
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0.nii.gz --mean --b0_thr $params.b0_thr_extract_b0
    """
}

dwi_for_bet
    .phase(b0_for_bet)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_b0_for_bet}

process bet_dwi {
    tag { "$sid" }
    cpus params.processes_bet_dwi

    input:
    set sid, file(dwi), file(b0) from dwi_b0_for_bet

    output:
    set sid, "${sid}__b0_bet_mask.nii.gz" into b0_mask_for_n4, b0_mask_for_crop
    set sid, "${sid}__dwi_bet.nii.gz" into dwi_for_n4
    set sid, "${sid}__b0_bet.nii.gz" into b0_for_n4, b0_for_crop

    script:
    dir_id = get_dir(sid)
    """
    antsBrainExtraction.sh -d 3 -a $b0 -e $template_dir_b0/b0_template.nii.gz\
        -o bet/ -m $template_dir_b0/b0_brain_probability_map.nii.gz\
        -f $template_dir_b0/b0_brain_registration_mask.nii.gz -k 1
    cp bet/BrainExtractionPriorWarped.nii.gz ${sid}__b0_bet_mask.nii.gz
    mrcalc $dwi ${sid}__b0_bet_mask.nii.gz -mult ${sid}__dwi_bet.nii.gz -quiet
    mrcalc $b0 ${sid}__b0_bet_mask.nii.gz -mult ${sid}__b0_bet.nii.gz -quiet
    """
}

dwi_for_n4
    .phase(b0_for_n4)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(b0_mask_for_n4)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_grad_b0_b0_mask_for_n4}

process n4_dwi {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(dwi), file(b0), file(b0_mask)\
        from dwi_grad_b0_b0_mask_for_n4

    output:
    set sid, "${sid}__dwi_n4.nii.gz" into dwi_for_crop

    script:
    dir_id = get_dir(sid)
    """
    N4BiasFieldCorrection -i $b0\
        -o [${sid}__b0_n4.nii.gz, bias_field_b0.nii.gz]\
        -c [300x150x75x50, 1e-6] -v 1
    scil_apply_bias_field_on_dwi.py $dwi bias_field_b0.nii.gz\
        ${sid}__dwi_n4.nii.gz --mask $b0_mask -f
    """
}

dwi_for_crop
    .phase(b0_mask_for_crop)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(b0_for_crop)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_and_b0_mask_b0_for_crop}

process crop_dwi {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(dwi), file(b0_mask), file(b0) from dwi_and_b0_mask_b0_for_crop

    output:
    set sid, "${sid}__dwi_crop.nii.gz", "${sid}__b0_mask_crop.nii.gz" into\
        dwi_and_b0_mask_for_denoise, dwi_mask_for_resample
    file "${sid}__b0_crop.nii.gz"

    script:
    dir_id = get_dir(sid)
    """
    scil_crop_volume.py $dwi ${sid}__dwi_crop.nii.gz -f\
        --output_bbox dwi_boundingBox.pkl -f
    scil_crop_volume.py $b0 ${sid}__b0_crop.nii.gz\
        --output_bbox b0_boundingBox.pkl -f
    scil_crop_volume.py $b0_mask ${sid}__b0_mask_crop.nii.gz\
        --input_bbox b0_boundingBox.pkl -f
    """
}

process resample_t1 {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(t1) from t1_for_resample

    output:
    set sid, "${sid}__t1_resample.nii.gz" into t1_for_bet

    script:
    dir_id = get_dir(sid)
    if(params.do_resample_t1)
        """
        scil_resample_volume.py $t1 ${sid}__t1_resample.nii.gz \
            --resolution $params.t1_resolution \
            --interp  $params.t1_interpolation
        """
    else
        """
        cp $t1 ${sid}__t1_resample.nii.gz
        """
}

process bet_t1 {
    tag { "$sid" }
    cpus params.processes_bet_t1

    input:
    set sid, file(t1) from t1_for_bet

    output:
    set sid, "${sid}__t1_bet.nii.gz" into t1_for_denoise
    set sid, "${sid}__t1_bet_mask.nii.gz" into t1_mask_for_denoise, t1_mask_for_crop

    script:
    dir_id = get_dir(sid)
    """
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$task.cpus
    antsBrainExtraction.sh -d 3 -a $t1 -e $template_dir_t1/t1_template.nii.gz\
        -o bet/ -m $template_dir_t1/t1_brain_probability_map.nii.gz \
        -f $template_dir_t1/t1_brain_registration_mask.nii.gz
    cp bet/BrainExtractionBrain.nii.gz ${sid}__t1_bet.nii.gz
    cp bet/BrainExtractionMask.nii.gz ${sid}__t1_bet_mask.nii.gz
    """
}

t1_for_denoise
    .phase(t1_mask_for_denoise)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{t1_and_mask_for_denoise}

process denoise_t1 {
    tag { "$sid" }
    cpus params.processes_denoise_t1

    input:
    set sid, file(t1), file(t1_mask) from t1_and_mask_for_denoise

    output:
    set sid, "${sid}__t1_denoised.nii.gz" into t1_for_n4

    script:
    dir_id = get_dir(sid)
    """
    scil_run_nlmeans.py $t1 ${sid}__t1_denoised.nii.gz 1 \
	    --mask $t1_mask --noise_est basic --processes $task.cpus -f
    """
}

process n4_t1 {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(t1) from t1_for_n4

    output:
    set sid, "${sid}__t1_n4.nii.gz" into t1_for_crop

    script:
    dir_id = get_dir(sid)
    """
    N4BiasFieldCorrection -i $t1\
        -o [${sid}__t1_n4.nii.gz, bias_field_t1.nii.gz]\
        -c [300x150x75x50, 1e-6] -v 1
    """
}

t1_for_crop
    .phase(t1_mask_for_crop)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{t1_and_mask_for_crop}

process crop_t1 {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(t1), file(t1_mask) from t1_and_mask_for_crop

    output:
    set sid, "${sid}__t1_bet_crop.nii.gz" into t1_for_reg
    set sid, "${sid}__t1_bet_mask_crop.nii.gz" into t1_mask_for_reg

    script:
    dir_id = get_dir(sid)
    """
    scil_crop_volume.py $t1 ${sid}__t1_bet_crop.nii.gz\
        --output_bbox t1_boundingBox.pkl -f
    scil_crop_volume.py $t1 ${sid}__t1_bet_mask_crop.nii.gz\
        --input_bbox t1_boundingBox.pkl -f
    """
}

process resample_dwi {
    tag { "$sid" }
    cpus 2

    input:
    set sid, file(dwi), file(mask) from dwi_mask_for_resample

    output:
    set sid, "${sid}__dwi_resample.nii.gz" into dwi_for_resample_b0,
        dwi_for_extract_dti_shell, dwi_for_extract_fodf_shell

    script:
    dir_id = get_dir(sid)
    if (params.do_resample_dwi)
        """
        scil_resample_volume.py $dwi \
            dwi_resample.nii.gz \
            --resolution $params.dwi_resolution \
            --interp  $params.dwi_interpolation
        scil_resample_volume.py $mask \
            mask_resample.nii.gz \
            --ref dwi_resample.nii.gz \
            --enforce_dimensions \
            --interp nn
        mrcalc dwi_resample.nii.gz mask_resample.nii.gz -mult ${sid}__dwi_resample.nii.gz -quiet
        """
    else
        """
        cp $dwi dwi_resample.nii.gz
        """
}

dwi_for_resample_b0
    .phase(bvals_for_resample_b0)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(bvecs_for_resample_b0)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_and_grad_for_resample_b0}

process resample_b0 {
    tag { "$sid" }
    cpus 2

    input:
    set sid, file(dwi), file(bval), file(bvec) from dwi_and_grad_for_resample_b0

    output:
    set sid, "${sid}__b0_resample.nii.gz" into b0_for_reg
    set sid, "${sid}__b0_mask_resample.nii.gz" into b0_mask_for_dti_metrics,
        b0_mask_for_fodf

    script:
    dir_id = get_dir(sid)
    """
    scil_extract_b0.py $dwi $bval $bvec ${sid}__b0_resample.nii.gz
    mrthreshold ${sid}__b0_resample.nii.gz ${sid}__b0_mask_resample.nii.gz\
        --abs 0.00001
    """
}

dwi_for_extract_dti_shell
    .phase(bvals_for_extract_dti_shell)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(bvecs_for_dti_shell)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_and_grad_for_extract_dti_shell}

process extract_dti_shell {
    tag { "$sid" }
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec) from\
        dwi_and_grad_for_extract_dti_shell

    output:
    set sid, "${sid}__dwi_dti.nii.gz", "${sid}__bval_dti",
        "${sid}__bvec_dti" into dwi_and_grad_for_dti_metrics

    script:
    dir_id = get_dir(sid)
    """
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.dti_shells ${sid}__dwi_dti.nii.gz \
        ${sid}__bval_dti ${sid}__bvec_dti -t $params.dwi_shell_tolerance -f
    """
}

dwi_and_grad_for_dti_metrics
    .phase(b0_mask_for_dti_metrics)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_and_grad_for_dti_metrics}

process dti_metrics {
    tag { "$sid" }
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask) from\
        dwi_and_grad_for_dti_metrics

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
    set sid, "${sid}__fa.nii.gz" into fa_for_reg

    script:
    dir_id = get_dir(sid)
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
    .phase(bvals_for_extract_fodf_shell)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(bvecs_for_fodf_shell)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{dwi_and_grad_for_extract_fodf_shell}

process extract_fodf_shell {
    tag { "$sid" }
    cpus 3

    input:
    set sid, file(dwi), file(bval), file(bvec) from\
        dwi_and_grad_for_extract_fodf_shell

    output:
    set sid, "${sid}__dwi_fodf.nii.gz", "${sid}__bval_fodf",
        "${sid}__bvec_fodf" into dwi_and_grad_for_fodf

    script:
    dir_id = get_dir(sid)
    """
    scil_extract_dwi_shell.py $dwi \
        $bval $bvec $params.fodf_shells ${sid}__dwi_fodf.nii.gz \
        ${sid}__bval_fodf ${sid}__bvec_fodf -t $params.dwi_shell_tolerance -f
    """
}

t1_for_reg
    .phase(t1_mask_for_reg)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(fa_for_reg)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(b0_for_reg)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{t1_fa_b0_for_reg}

process register_t1 {
    tag { "$sid" }
    cpus params.processes_registration

    input:
    set sid, file(t1), file(t1_mask), file(fa), file(b0) from t1_fa_b0_for_reg

    output:
    set sid, "${sid}__t1_warp.nii.gz" into t1_for_seg
    file "${sid}__output0GenericAffine.mat"
    file "${sid}__output1InverseWarp.nii.gz"
    file "${sid}__output1Warp.nii.gz"
    file "${sid}__t1_mask_warp.nii.gz"

    script:
    dir_id = get_dir(sid)
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
    cp outputWarped.nii.gz ${sid}__t1_warp.nii.gz
    cp output0GenericAffine.mat ${sid}__output0GenericAffine.mat
    cp output1InverseWarp.nii.gz ${sid}__output1InverseWarp.nii.gz
    cp output1Warp.nii.gz ${sid}__output1Warp.nii.gz
    antsApplyTransforms -d 3 -i $t1_mask -r ${sid}__t1_warp.nii.gz \
        -o ${sid}__t1_mask_warp.nii.gz -n NearestNeighbor \
        -t ${sid}__output1Warp.nii.gz ${sid}__output0GenericAffine.mat
    """
}

process segment_tissues {
    tag { "$sid" }
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
    dir_id = get_dir(sid)
    """
    fast -t $params.type_of_image -n $params.number_of_tissue\
         -H $params.spatial_smoothness\
         -I $params.number_of_iter_for_bias_field\
         -l $params.lowpass -g -o t1.nii.gz $t1
    cp t1_seg_2.nii.gz ${sid}__mask_wm.nii.gz
    cp t1_seg_1.nii.gz ${sid}__mask_gm.nii.gz
    cp t1_seg_0.nii.gz ${sid}__mask_csf.nii.gz
    cp t1_pve_2.nii.gz ${sid}__map_wm.nii.gz
    cp t1_pve_1.nii.gz ${sid}__map_gm.nii.gz
    cp t1_pve_0.nii.gz ${sid}__map_csf.nii.gz
    """
}

dwi_and_grad_for_fodf
    .phase(b0_mask_for_fodf)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .phase(fa_md_for_fodf)
    .map{ch1, ch2 -> [*ch1, ch2[1], ch2[2]] }
    .set{dwi_b0_metrics_for_fodf}

process fodf_metrics {
    tag { "$sid" }
    cpus params.processes_fodf

    input:
    set sid, file(dwi), file(bval), file(bvec), file(b0_mask), file(fa),
        file(md) from dwi_b0_metrics_for_fodf

    output:
    set sid, "${sid}__fodf.nii.gz" into fodf_for_tracking
    file "${sid}__peaks.nii.gz"
    file "${sid}__peak_indices.nii.gz"
    file "${sid}__afd_max.nii.gz"
    file "${sid}__afd_total.nii.gz"
    file "${sid}__afd_sum.nii.gz"
    file "${sid}__nufo.nii.gz"

    script:
    dir_id = get_dir(sid)
    """
    scil_compute_ssst_frf.py $dwi $bval $bvec ${sid}__frf.txt --mask $b0_mask\
        --fa $params.fa --min_fa $params.min_fa --min_nvox $params.min_nvox\
        --roi_radius $params.roi_radius
        
    scil_compute_fodf.py $dwi $bval $bvec ${sid}__frf.txt --sh_order $params.sh_order\
        --basis $params.basis --b0_threshold $params.b0_thr_extract_b0 \
        --mask $b0_mask --fodf ${sid}__fodf.nii.gz --peaks ${sid}__peaks.nii.gz\
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

process pft_maps {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(wm), file(gm), file(csf) from map_wm_gm_csf_for_pft_maps

    output:
    set sid, "${sid}__map_include.nii.gz", "${sid}__map_exclude.nii.gz" into pft_maps_for_tracking
    set sid, "${sid}__interface.nii.gz" into interface_for_seeding_mask

    script:
    dir_id = get_dir(sid)
    """
    scil_compute_maps_for_particle_filter_tracking.py $wm $gm $csf \
    --include ${sid}__map_include.nii.gz --exclude ${sid}__map_exclude.nii.gz\
        --interface ${sid}__interface.nii.gz -f
    """
}


wm_mask_for_seeding_mask
    .phase(interface_for_seeding_mask)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{wm_interface_for_seeding_mask}
process seeding_mask {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(wm), file(interface_mask) from wm_interface_for_seeding_mask

    output:
    set sid, "${sid}__seeding_mask.nii.gz" into seeding_mask_for_tracking

    script:
    dir_id = get_dir(sid)
    if (params.wm_seeding)
        """
        scil_mask_math.py union $wm $interface_mask ${sid}__seeding_mask.nii.gz
        """
    else
        """
        cp $interface_mask ${sid}__seeding_mask.nii.gz
        """
}

fodf_for_tracking
    .phase(pft_maps_for_tracking)
    .map{ch1, ch2 -> [*ch1, ch2[1], ch2[2]] }
    .phase(seeding_mask_for_tracking)
    .map{ch1, ch2 -> [*ch1, ch2[1]] }
    .set{fodf_maps_for_tracking}

process tracking {
    tag { "$sid" }
    cpus params.processes_tracking

    input:
    set sid, file(fodf), file(include), file(exclude), file(seed) from\
        fodf_maps_for_tracking

    output:
    file "${sid}__tracking.trk"

    script:
    dir_id = get_dir(sid)
    if (params.do_compress_streamlines)
        """
        scil_compute_particle_filter_tracking.py $fodf $seed\
            $include $exclude ${sid}__tracking.trk --algo $params.algo\
            --$params.seeding $params.nbr_seeds --random $params.random --step $params.step\
            --rk_order $params.rk_order --theta $params.theta --maxL_no_dir $params.maxL_no_dir\
            --sfthres $params.sfthres --sfthres_init $params.sfthres_init --minL $params.minL\
            --maxL $params.maxL --sh_interp $params.sh_interp --mask_interp $params.mask_interp\
            --particles $params.particles --back $params.back --front $params.front\
            --pft_theta $params.pft_theta --processes $task.cpus --compress $params.compress_value
        """
    else
        """
        scil_compute_particle_filter_tracking.py $fodf $seed\
            $include $exclude ${sid}__tracking.trk --algo $params.algo\
            --$params.seeding $params.nbr_seeds --random $params.random --step $params.step\
            --rk_order $params.rk_order --theta $params.theta --maxL_no_dir $params.maxL_no_dir\
            --sfthres $params.sfthres --sfthres_init $params.sfthres_init --minL $params.minL\
            --maxL $params.maxL --sh_interp $params.sh_interp --mask_interp $params.mask_interp\
            --particles $params.particles --back $params.back --front $params.front\
            --pft_theta $params.pft_theta --processes $task.cpus
        """
}

def get_dir(sid) {
    String my_new_str = sid.replaceAll("_-_", "/")
    return my_new_str
}