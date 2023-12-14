import SimpleITK as sitk
import numpy as np
from radiomics.imageoperations import checkMask, getMask
import os

for mask_folder in snakemake.input.mask_folders:
    maskname = os.path.join(mask_folder, '1-1.dcm')
    mask = sitk.ReadImage(maskname)
    mask = sitk.Cast(mask, sitk.sitkInt32)
    mask = sitk.DICOMOrient(mask, 'LPS')
    name = mask_folder+'.nrrd'
    sitk.WriteImage(mask, name)

    