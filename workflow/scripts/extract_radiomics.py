import csv
import SimpleITK as sitk
from radiomics import featureextractor

extractor = featureextractor.RadiomicsFeatureExtractor()
extractor.loadParams(
    '/mnt/home/fesslerl/890-402/cmse_890_602_fessler_final/workflow/env/radiomics_settings.yml')

print('Extraction parameters:\n\t', extractor.settings)
print('Enabled filters:\n\t', extractor.enabledImagetypes)
print('Enabled features:\n\t', extractor.enabledFeatures)

for (cropped_mask_filename,
     cropped_scan_filename,
     radiomics_filename) in zip(
         snakemake.input.cropped_mask_filenames,
         snakemake.input.cropped_scan_filenames,
         snakemake.output.radiomics_filenames):

    mask = sitk.ReadImage(cropped_mask_filename)
    scan = sitk.ReadImage(cropped_scan_filename)

    # Remember that masks are provided as inverse masks where 0 indicates tumor masking. We will need to invert this mask ourselves.
    mask_inverted = (mask == 0)

    features = extractor.execute(scan, mask_inverted)

    # features is an ordered dictionary. We will write this output to a csv,
    # being careful to preserve order.
    with open(radiomics_filename, 'a') as outputFile:
        writer = csv.writer(outputFile, lineterminator='\n')
        headers = list(features.keys())
        writer.writerow(headers)

        row = []
        for h in headers:
            row.append(features.get(h, "N/A"))
        writer.writerow(row)
