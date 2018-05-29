Human pipeline
==============

Do all the processing. The different steps are explained below. The input subjects
must all have the same encoding scheme (gradient table).

WARNING: The results could take 4 times more disk space than the input.

If you use this pipeline, please cite:

```
Kurtzer GM, Sochat V, Bauer MW (2017)
Singularity: Scientific containers for mobility of compute.
PLoS ONE 12(5): e0177459. https://doi.org/10.1371/journal.pone.0177459

P. Di Tommaso, et al. Nextflow enables reproducible computational workflows.
Nature Biotechnology 35, 316–319 (2017) doi:10.1038/nbt.3820
```

Steps
-----

######Diffusion processes
- preliminary DWI brain extraction (ANTs)
- denoise DWI (Mrtrix3)
- topup (FSL)
- eddy (FSL)
- extract B0 (Scilpy)
- dwi brain extraction (ANTs)
- N4 DWI (ANTs)
- crop DWI (Scilpy)
- upsample DWI (Scilpy)
- upsample B0 (Scilpy)

######T1 processes
- denoise T1 (Scilpy)
- N4 T1 (ANTs)
- resample T1 (Scilpy)
- T1 brain extraction (ANTs)
- crop T1 (Scilpy)
- registration on diffusion (ANTs)
- tissue segmentation (FSL)

######Metrics processes
- extract DTI shells (configurable in `nextflow.config`) (Scilpy)
- dti metrics (Scilpy)
- extract fODF shells (configurable in `nextflow.config`) (Scilpy)
- compute fiber response function (Scilpy)
- mean fiber response function (Scilpy)
- fODF and fODF-derived metrics (Scilpy)
- particle filter tractography maps (Scilpy)
- seeding mask (Scilpy)
- particle filter tracking (Scilpy)

Singularity
-----------
The singularity is available on
[singularity-human](https://bitbucket.org/sciludes/singularity-human/src/master/) repository
 
Requirements
------------

- [Nextflow](https://www.nextflow.io)

The file _singularity_human.def_ contains most of the requirements commands.

- Scilpy
- Dipy
- Mrtrix3
- ANTS
- FSL with the eddy_openmp
[patch](https://fsl.fmrib.ox.ac.uk/fsldownloads/patches/eddy-patch-fsl-5.0.9/centos6/)

Notes
-----

The _scilpy/scripts_ folder should be in your PATH environment variable if
Singularity is not used.

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`