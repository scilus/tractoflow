TractoFlow pipeline
===================

The TractoFlow pipeline is a fully automated and reproducible dMRI processing pipeline.
TractoFlow takes raw DWI, b-values, b-vectors, T1 weighted image (and a reversed
phase encoded b=0 if available) to process DTI, fODF metrics and a whole brain tractogram.

Documentation:
--------------

TractoFlow documentation is available here:

[https://tractoflow-documentation.readthedocs.io](https://tractoflow-documentation.readthedocs.io)

This documentation presents how to install and launch TractoFlow on a local computer and a High Performance Computer.

If you are a user and NOT A DEVELOPER, we HIGHLY RECOMMEND following the instructions on the documentation website.


Support:
--------

For questions or help on TractoFlow, please create a topic on Neurostar with the `tractoflow` tag [https://neurostars.org/tags/tractoflow](https://neurostars.org/tags/tractoflow)

Before posting, please check the topics already created. [https://neurostars.org/tags/tractoflow](https://neurostars.org/tags/tractoflow)


Singularity
-----------
If you are on Linux, we recommend using the Singularity container to run TractoFlow

Prebuild Singularity images are available here:

[http://scil.usherbrooke.ca/en/containers_list/](http://scil.usherbrooke.ca/en/containers_list/)

FOR DEVELOPERS: The Singularity repository is available here:
[singularity-TractoFlow](https://github.com/scilus/singularity-tractoflow)

Docker
------
If you are on MacOS or Windows, we recommend using the Docker container to run TractoFlow

Prebuilt Docker images are available here:

[http://scil.usherbrooke.ca/en/containers_list/](http://scil.usherbrooke.ca/en/containers_list/)

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

To refer to TractoFlow, please see `References` section on [TractoFlow documentation](https://tractoflow-documentation.readthedocs.io/en/latest/reference/references.html)
