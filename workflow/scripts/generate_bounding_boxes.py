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
