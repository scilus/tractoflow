Human pipeline
--------------------------

These pipelines process single-shell and multi-shell datasets. There are 2 pipelines:

- _preprocess_stride_flip/main.nf_ : Modify strides and flip gradients. This step is not
mandatory if the strides and flips are already correct. 

- _pipeline_multi_shell/main.nf_ : Do all the processing. The different steps are
explained below. The input dataset must have
the same gradient scheme.

Requirements
------------

- [Nextflow](https://www.nextflow.io)
- Scilpy
- Dipy
- Mrtrix3
- ANTS
- FSL with the eddy_openmp [patch](https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/)

The file _singularity/Singularity.def_ contains most of the requirements commands.