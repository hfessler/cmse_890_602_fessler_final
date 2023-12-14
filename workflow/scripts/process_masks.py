import os
import SimpleITK as sitk

for mask_folder in snakemake.input.mask_folders:

    # ISPY2 masks are segmentation dicom objects,
    # and are only 1 .dcm file long.
    mask_filename = os.path.join(mask_folder, '1-1.dcm')
    mask = sitk.ReadImage(mask_filename)

    # Masks should be cast to 32 bit integers and put in standard
    # orientation (right to Left, anterior to Posterior, inferior to Superior)
    mask = sitk.Cast(mask, sitk.sitkInt32)
    mask = sitk.DICOMOrient(mask, 'LPS')

    new_filename = mask_folder+'.nrrd'
    sitk.WriteImage(mask, new_filename)
