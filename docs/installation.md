# Minimal Package Requirements

To run this workflow, you will need **snakemake**. This workflow was built with snakemake **7.32.4**.

To test this workflow, you will need **pytest**. This workflow was tested with pytest ***7.4.3***. 

I have not had the opportunity to check if older versions of either software cause errors. 


# Using your own environment

All further dependencies are handled via a conda environment created by the workflow itself. If you would instead like to use your own environment, the key dependencies are:

```
python >= 3.9
tcia-utils >= 1.8
pyradiomics >= 3.0
simpleitk >= 2.3
pandas >= 2.1
```

Package versions below these values have not been tested and may or may not function. 

Be aware that without access to the conda-forge channel that building an environment with snakemake via conda sometimes builds with python 3.6. tcia-utils only has support for python 3.7 or greater. For building your own environment, you will have to be careful that this requirement is met. 

