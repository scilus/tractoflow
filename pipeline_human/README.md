Human pipeline
=================

Do all the processing. The different steps are
explained below. The input subjects must all have the same encoding scheme (gradient table).

WARNING: The results could take 4 times more disk space than the input. 

Steps
------------

######Diffusion processes
- preliminary dwi brain extraction (ANTs)
- denoise dwi (Mrtrix3)
- topup (FSL)
- eddy (FSL)
- extract b0 (Scilpy)
- dwi brain extraction (ANTs)
- N4 dwi (ANTs)
- crop dwi (Scilpy)
- upsample dwi (Scilpy)
- upsample b0 (Scilpy)

######T1 processes
- resample t1 (Scilpy)
- t1 brain extraction (ANTs)
- N4 t1 (ANTs)
- crop t1 (Scilpy)
- denoise t1 (Scilpy)
- registration on diffusion (ANTs)
- tissue segmentation (FSL)

######Metrics processes
- extract dti shells (b-values < 1200) (Scilpy)
- dti metrics (Scilpy)
- extract fodf shells (b-values > 700) (Scilpy)
- compute fiber response function (Scilpy)
- mean fiber response function (Scilpy)
- fodf metrics (Scilpy)
- particle filter tractography maps (Scilpy)
- seeding_mask (Scilpy)
- particle filter tracking (Scilpy)

Singularity
-----
The image for singularity can be built using _singularity/Singularity.def_ with the command:
```singularity build image_name.img Singularity.def```. It could be used to run
the human pipeline with the option ```-with-singularity image_name.img```
 of nextflow.
 To build the singularity, please run the command from the directory ```singularity/```. Otherwise, ```scilpy.tar``` is not found.
 
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