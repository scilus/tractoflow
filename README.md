TractoFlow pipeline
===================

TractoFlow pipeline is a fully automated and reproducible dMRI pipeline.
TractoFlow takes raw DWI, b-values, b-vectors, T1 weighted image (and a reversed
phase encoded b=0 if available) to process DTI, fODF metrics and a whole brain tractogram.

Documentation:
--------------

TractoFlow documentation is available here:

[https://tractoflow-documentation.readthedocs.io](https://tractoflow-documentation.readthedocs.io)

This documentation present who to install and launch TractoFlow on a local computer and an High Performance Computer.

If you are a user and NOT A DEVELOPPER, we HIGHLY RECOMMEND to follow the instructions on the documentation website.

Singularity
-----------
If you are on Linux, please use the Singularity container the run TractoFlow

Please download Singularity container built here:

[http://scil.dinf.usherbrooke.ca/en/containers_list/](http://scil.dinf.usherbrooke.ca/en/containers_list/)

FOR DEVELOPPERS: The Singularity repository is available here:
[singularity-TractoFlow](https://github.com/scilus/singularity-tractoflow)

Docker
------
If you are on MacOS or Windows, please use the Docker container the run TractoFlow

Please download Docker container built here:

[http://scil.dinf.usherbrooke.ca/en/containers_list/](http://scil.dinf.usherbrooke.ca/en/containers_list/)

Steps
-----

Diffusion processes
- preliminary DWI brain extraction (FSL)
- denoise DWI (Mrtrix3)
- topup (FSL)
- eddy (FSL)
- extract B0 (Scilpy)
- dwi brain extraction (FSL)
- N4 DWI (ANTs)
- crop DWI (Scilpy)
- upsample DWI (Scilpy)
- upsample B0 (Scilpy)

T1 processes
- denoise T1 (Scilpy)
- N4 T1 (ANTs)
- resample T1 (Scilpy)
- T1 brain extraction (ANTs)
- crop T1 (Scilpy)
- registration on diffusion (ANTs)
- tissue segmentation (FSL)

Metrics processes
- extract DTI shells (configurable in `nextflow.config`) (Scilpy)
- dti metrics (Scilpy)
- extract fODF shells (configurable in `nextflow.config`) (Scilpy)
- compute fiber response function (Scilpy)
- mean fiber response function (Scilpy)
- fODF and fODF-derived metrics (Scilpy)
- particle filter tractography maps (Scilpy)
- seeding mask (Scilpy)
- particle filter tracking (Scilpy)

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`

References
----------

To refer to TratoFlow, please see `References` section on [TractoFlow documentation](https://tractoflow-documentation.readthedocs.io)
