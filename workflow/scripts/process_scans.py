import SimpleITK as sitk

for scan_folder in snakemake.input.scan_folders:

    # ISPY2 scan images have a variable number of dcm files,
    # one for each slice of the image
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(scan_folder)
    reader.SetFileNames(dicom_names)
    scan = reader.Execute()

    # Scans should be put in standard orientation (right to Left,
    # anterior to Posterior, inferior to Superior)
    scan = sitk.DICOMOrient(scan, 'LPS')

    new_filename = scan_folder+'.nrrd'
    sitk.WriteImage(scan, new_filename)
