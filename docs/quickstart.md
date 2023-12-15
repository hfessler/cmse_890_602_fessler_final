# Installing snakemake

You will need snakemake to run this workflow. To create an environment with conda, you can run the following command:

```bash
conda create -n radiomics_of_ISPY -c conda-forge bioconda::snakemake=7.32.4
```

Be sure to activate this enviroment before proceeding. 

```bash
conda activate radiomics_of_ISPY
```

# Running the workflow
From here you can run the workflow. 

```bash
snakemake --use-conda --cores 2
```
The use-conda flag will have the workflow create its own conda environment. You are free to change the number of cores requested. 

# Changing global variables
The key variables to change for the workflow are the collection name, the scan name, the mask name, and patient id(s), the times, and the permanent storage location. All of these are global variables at the very top of the Snakefile. 

For example, if I wanted to collect radiomics from two patients at first and last observation from the sister collection ACRIN 6698 who were uni-laterally cropped and store these in my scratch space, I would change the header of the Snakefile to 
```
# Collection and scan and mask identifying names.
COLLECTION = 'ACRIN-6698'
SCAN_NAME = 'ISPY2: VOLSER: uni-lateral cropped: SER'
MASK_NAME = 'ISPY2: VOLSER: uni-lateral cropped: Analysis Mask'

# Images in ISPY2 are unique for a specific patient and a specific time, because each patient is imaged multiple times. Times in ISPY2 are a simple ordering of images, not a specific calendar time.
PATIENT_IDS = ['ACRIN-6698-102212', 'ACRIN-6698-103939']
TIMES = [0,3]

# Path to a permanent location to store the outputs of this workflow.
PERM_STORAGE='/mnt/scratch/username'
```

