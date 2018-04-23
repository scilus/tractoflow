Human pipeline
=================

Do all the processing. The different steps are
explained below. The input dataset must have
the same gradient scheme.

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
- upsample dwi
- upsample b0

######T1 processes
- resample t1
- bet t1
- n4 t1
- crop t1
- denoise t1
- registration on diffusion
- tissue segmentation

######Metrics processes
- extract dti shell
- dti metrics
- extract fodf shell
- compute frf
- mean frf
- fodf metrics
- pft maps
- seeding_mask
- tracking

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

Usage
-----

See *USAGE* or run ```nextflow run main.nf --help```