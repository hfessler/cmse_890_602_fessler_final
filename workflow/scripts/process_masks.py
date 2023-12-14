import SimpleITK as sitk
import numpy as np
from radiomics.imageoperations import checkMask, getMask
import os

for i, mask_folder in enumerate(snakemake.input.mask_folders):
    maskname = os.path.join(mask_folder, '1-1.dcm')
    mask_sitk = sitk.ReadImage(maskname)
    mask_sitk = sitk.Cast(mask_sitk, sitk.sitkInt32)
    # ISPY masks are inverse masks, so tumor is marked by label 0
    mask_inverted = (mask_sitk==0)
    scan = sitk.ReadImage(snakemake.input.scans[i])

    mask_inverted = sitk.DICOMOrient(mask_inverted, 'LPS')
    scan = sitk.DICOMOrient(scan, 'LPS')

    bounding_box, _ = checkMask(scan, mask_inverted)

    sitk.WriteImage(mask_sitk, snakemake.output.masks[i])    
    np.save(snakemake.output.bounding_boxes[i], bounding_box)