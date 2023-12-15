"""Generate a bounding box for a ROI

    Generate a bounding box for a mask scan pair. 

    Args:
        scan_filenames: List of strings of filenames of scan files. Scan files are expected to be in NRRD format. 
        mask_filenames: List of strings of filenames of mask files. Mask files are expected to be in NRRD format. 
        bounding_box_filename: List of strings of filenames for bounding boxes to be saved as.

    Extract the bounding box capturing the region of interest (ROI) defined by an inverse 
    mask and save the output in a numpy array. 

    Note that the lists are expected to be ordered similarly. That is, the nth element 
    of each list corresponds to the same imaging event. Failure to do so will result in 
    mislabeling or comparison of unrelated images.  

    Further these masks are expected to be inverse masks, where tumor regions are marked 
    with value 0. This is in contrast to the more common masks where tumor regions are 
    marked 1. 
    """


import numpy as np
import SimpleITK as sitk
from radiomics.imageoperations import checkMask

for (scan_filename,
     mask_filename,
     bounding_box_filename) in zip(
         snakemake.input.scan_filenames,
         snakemake.input.mask_filenames,
         snakemake.output.bounding_box_filenames):

    scan = sitk.ReadImage(scan_filename)
    mask = sitk.ReadImage(mask_filename)

    bounding_box, _ = checkMask(scan, mask, label=0)

    np.save(bounding_box_filename, bounding_box)
