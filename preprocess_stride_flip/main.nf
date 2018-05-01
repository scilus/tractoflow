#!/usr/bin/env nextflow

params.root = false
params.subject = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    cpu_count = Runtime.runtime.availableProcessors()
    bindings = ["flipx":"$params.flipx", "flipy":"$params.flipy",
                "flipz":"$params.flipz"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Preprocess pipeline"
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
        .fromFilePairs("$root/**/{bval,bvec,dwi.nii.gz,t1.nii.gz}",
                       size: 4,
                       maxDepth:2,
                       flat: true) { it.parent.parent.name + "_-_" + it.parent.name}
    Channel
    .fromPath("$root/*rev_b0.nii.gz",
                    maxDepth:2)
    .map{ch1 -> [ch1.parent.parent.name + "_-_" + ch1.parent.name, ch1]}
    .into{rev_b0; check_rev_b0}
    }
else if (params.subject){
    log.info "Input: $params.subject"
    subject = file(params.subject)
    in_data = Channel
        .fromFilePairs("$subject/{bval,bvec,dwi.nii.gz,t1.nii.gz}",
                       size: 4,
                       maxDepth:2,
                       flat: true) { it.parent.name}
    Channel
    .fromPath("$subject/*rev_b0.nii.gz",
                    maxDepth:1)
    .map{ch1 -> [ch1.parent.name, ch1]}
    .into{rev_b0; check_rev_b0}
    }

(dwi_t1_for_convert, gradients_for_flip) = in_data
    .map{sid, bvals, bvecs, dwi, t1 -> [tuple(sid, dwi, t1),
                                        tuple(sid, bvals, bvecs)]}
    .separate(2)

check_rev_b0.count().set{ rev_b0_counter }

process correct_stride {
    tag { "$sid" }
    cpus 2

    input:
    set sid, file(dwi), file(t1) from dwi_t1_for_convert

    output:
    file "dwi.nii.gz"
    file "t1.nii.gz"

    script:
    dir_id = get_dir(sid)
    """
    mrconvert $dwi dwi.nii.gz -stride 1,2,3,4 -force -quiet
    mrconvert $t1 t1.nii.gz -stride 1,2,3,4 -force -quiet
    """
}

process correct_stride_rev_b0 {
    tag { "$sid" }
    cpus 1

    input:
    set sid, file(rev) from rev_b0

    output:
    file "rev_b0.nii.gz"

    when:
    check_rev_b0 == 1

    script:
    dir_id = get_dir(sid)
    """
    mrconvert $rev rev_b0.nii.gz -stride 1,2,3,4 -force -quiet
    """
}

process flip_gradients {
    tag { "$sid" }
    cpus 1
    echo true

    input:
    set sid, file(bval_in), file(bvec_in) from gradients_for_flip

    output:
    file "bval"
    file "bvec"

    script:
    dir_id = get_dir(sid)
    """
    touch bval
    """

    String flip = ""
    if (params.flipx)
        flip += " x"
    if (params.flipy)
        flip += " y"
    if (params.flipz)
        flip += " z"
    if (params.flipx || params.flipy || params.flipz)
        """
        scil_flip_grad.py $bvec_in bvec $flip --fsl -f
        """
    else
        """
        touch bvec
        """
}

def get_dir(sid) {
    String my_new_str = sid.replaceAll("_-_", "/")
    return my_new_str
}