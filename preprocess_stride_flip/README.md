Penthera pipeline
=================

These pipelines process datasets similar to Penthera 3T. There are 2 pipelines:

- _preprocess_penthera/main.nf_ : Modify strides and flip gradients. This step is not
mandatory if the strides and flips are already correct. The input dataset must have
the same gradient scheme.

- _pipeline_penthera/main.nf_ : Do all the processing. The different steps are
explained below.

Steps
------------

######Diffusion processes
- preliminary bet
- denoise dwi
- eddy
- extract b0
- bet dwi
- n4 dwi
- crop dwi
- upsample dwi (1mm iso)
- upsample b0 (1mm iso)

######T1 processes
- resample t1 (1mm iso)
- bet t1
- n4 t1
- crop t1
- denoise t1
- registration on diffusion
- tissue segmentation

######Metrics processes
- extract dti shell (b=0, 300, 1000)
- dti metrics
- extract fodf shell (b=0, 2000)
- fodf metrics
- pft maps
- tracking (interface seeding)

Singularity
-----
The image for singularity can be built using _singularity/Singularity.def_ with the command:
```singularity build image_name.img Singularity.def```. It could be used to run
the penthera pipelines with the option ```-with-singularity image_name.img```
 of nextflow.
 
 
Requirements
------------

- [Nextflow](https://www.nextflow.io)
- Scilpy
- Dipy
- Mrtrix3
- ANTS
- FSL with the eddy_openmp [patch](https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/)

The file _singularity/Singularity.def_ contains most of the requirements commands.

Notes
-----

The _scilpy/scripts_ folder should be in your PATH environment variable.


The _penthera-nf/scripts_ folder should be in your PATH environment variable.

Usage
-----

See *USAGE* or run ```nextflow run main.nf --help```