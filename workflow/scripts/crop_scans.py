import SimpleITK as sitk
import numpy as np
from radiomics.imageoperations import cropToTumorMask

for mask_filename, scan_filename, bounding_box_filename, cropped_scan_filename, cropped_mask_filename in zip(snakemake.input.masks, snakemake.input.scans, snakemake.input.bounding_boxes, snakemake.output.cropped_scans, snakemake.output.cropped_masks):
    mask = sitk.ReadImage(mask_filename)
    scan = sitk.ReadImage(scan_filename)
    bounding_box = np.load(bounding_box_filename)
    cropped_scan, cropped_mask = cropToTumorMask(scan, mask, bounding_box)
    sitk.WriteImage(cropped_scan, cropped_scan_filename)
    sitk.WriteImage(cropped_mask, cropped_mask_filename)
