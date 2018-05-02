Human pipeline
--------------------------

The image for singularity can be built using _singularity/Singularity.def_ with the command:
```singularity build image_name.img Singularity.def```. It can be used to run
the human pipeline with the option ```-with-singularity image_name.img```
 of nextflow.
 To build the singularity, please run the command from the directory ```singularity/```. Otherwise, ```scilpy.tar``` is not found.
 