import SimpleITK as sitk  

for scan_folder in snakemake.input:
    print(scan_folder)
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(scan_folder)
    reader.SetFileNames(dicom_names)
    image = reader.Execute()
    
    name = scan_folder+'.nrrd'
    sitk.WriteImage(image, name)