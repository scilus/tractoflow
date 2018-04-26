Human pipeline
=================

Modify strides and flip gradients. This step is not
mandatory if the strides of the DWI is [1 2 3 4].

Singularity
-----
The image for singularity can be built using _singularity/Singularity.def_ with the command:
```singularity build image_name.img Singularity.def```. It could be used to run
the preprocess pipeline with the option ```-with-singularity image_name.img```
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