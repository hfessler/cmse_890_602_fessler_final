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
