"""Extract radiomics from a scan and mask pair

    Extract radiomics from a paired scan and inverse mask of the same size and orientation. 

    Args:
        cropped_mask_filenames: List of strings of filenames of cropped mask files. Mask 
        files are expected to be in NRRD format. 
        cropped_scan_filenames: List of strings of filenames of cropped scan files. Scan 
        files are expected to be in NRRD format. 
        radiomics_filenames: List of strings of filenames for extracted radiomics from mask scan pairs. 

    Extracts radiomics for each scan mask pair and saves the results to a csv.

    Note that the lists are expected to be ordered similarly. That is, the nth element 
    of each list corresponds to the same imaging event. Failure to do so will result in 
    mislabeling or comparison of unrelated images.  

    Further these masks are expected to be inverse masks, were tumor regions are marked 
    with value 0. This is in contrast to the more common masks where tumor regions are 
    marked 1. 
    """


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

    # Remember that masks are provided as inverse masks where 0
    # indicates tumor masking. We will need to invert this mask ourselves.
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
