import SimpleITK as sitk
import csv
from radiomics import featureextractor, getFeatureClasses

extractor = featureextractor.RadiomicsFeatureExtractor()
extractor.loadParams('/mnt/home/fesslerl/890-402/cmse_890_602_fessler_final/workflow/env/radiomics_settings.yml')

print('Extraction parameters:\n\t', extractor.settings)
print('Enabled filters:\n\t', extractor.enabledImagetypes)
print('Enabled features:\n\t', extractor.enabledFeatures)

for cropped_mask_filename, cropped_scan_filename, radiomics_filename in zip(snakemake.input.cropped_masks, snakemake.input.cropped_scans, snakemake.output.radiomics):

    mask = sitk.ReadImage(cropped_mask_filename)
    scan = sitk.ReadImage(cropped_scan_filename)
    # Remember that masks are provided as inverse masks where 0 indicates tumor masking
    mask_inverted = (mask==0)

    features = extractor.execute(scan, mask)
    with open(radiomics_filename, 'a') as outputFile:
        writer = csv.writer(outputFile, lineterminator='\n')
        headers = list(features.keys())
        writer.writerow(headers)
        row = []
        for h in headers:
            row.append(features.get(h, "N/A"))
        writer.writerow(row)