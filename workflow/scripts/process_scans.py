import SimpleITK as sitk

for scan_folder in snakemake.input.scan_folders:
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(scan_folder)
    reader.SetFileNames(dicom_names)
    scan = reader.Execute()
    scan = sitk.DICOMOrient(scan, 'LPS')

    name = scan_folder+'.nrrd'
    sitk.WriteImage(scan, name)
