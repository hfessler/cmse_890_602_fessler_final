"""Crop scans to area of interest

    Takes in a scan and mask and crops the scan and mask to the bounding box. 

    Args:
        mask_filenames: List of strings of filenames of mask files. Mask files are expected to be in NRRD format. 
        scan_filenames:  List of strings of filenames of scan files. Scan files are expected to be in NRRD format and must match shape and orientation of corresponding mask.  
        bounding_box_filenames: List of strings of filenames of bounding boxes. Expected to be in npy format. Should be of the format returned by radiomics.imageoperations.CheckMask. 
        cropped_scan_filenames: List of strings of filenames to save the cropped scans under. 
        cropped_mask_filenames: List of strings of filenames to save the cropped masks under. 

    
    Saves the cropped scans and cropped masks under the provided filenames. 

    Note that the lists are expected to be ordered similarly. That is, the nth element 
    of each list corresponds to the same imaging event. Failure to do so will result in 
    mislabeling or comparison of unrelated images.  
    The cropped masks are marked as temporary in the workflow. They will not be placed into 
    permanent storage, and will be removed when not needed. 
    """


import numpy as np
import SimpleITK as sitk
from radiomics.imageoperations import cropToTumorMask

for (mask_filename,
     scan_filename,
     bounding_box_filename,
     cropped_scan_filename,
     cropped_mask_filename) in zip(
        snakemake.input.mask_filenames,
        snakemake.input.scan_filenames,
        snakemake.input.bounding_box_filenames,
        snakemake.output.cropped_scan_filenames,
        snakemake.output.cropped_mask_filenames):

    mask = sitk.ReadImage(mask_filename)
    scan = sitk.ReadImage(scan_filename)
    bounding_box = np.load(bounding_box_filename)

    cropped_scan, cropped_mask = cropToTumorMask(scan, mask, bounding_box)

    sitk.WriteImage(cropped_scan, cropped_scan_filename)
    sitk.WriteImage(cropped_mask, cropped_mask_filename)
