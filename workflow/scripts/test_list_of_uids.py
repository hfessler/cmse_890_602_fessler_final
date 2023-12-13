import pandas as pd
from tcia_utils import nbia

uids=pd.read_csv(snakemake.input.uids)

report = open(snakemake.output.uid_report, "a")

report.write(f"{len(uids)} studies found with {snakemake.params.scan} series from {snakemake.params.collection}.\n")

report.write(f"{uids['mask_UID'].isna().sum()} ({round(uids['mask_UID'].isna().sum()/len(uids)*100, 2)}%) studies do not have an matching {snakemake.params.mask} series.\n")

uids_with_mask = uids.dropna()
uids_with_mask = uids_with_mask.reset_index() 

mismatched_scan_mask = []
for index, row in uids_with_mask.iterrows():
    scan_UID = row['scan_UID']
    mask_UID = row['mask_UID']

    scan_metadata = nbia.getSeriesMetadata(scan_UID)[0]
    mask_metadata = nbia.getSeriesMetadata(mask_UID)[0]

    # studies in ISPY are unique for subject_ID and timing of study
    scan_time = scan_metadata['Study Description'][-1]
    scan_subject_id = scan_metadata['Subject ID'][-6:]

    mask_time = mask_metadata['Study Description'][-1]
    mask_subject_id = mask_metadata['Subject ID'][-6:]

    if mask_time != scan_time or scan_subject_id != mask_subject_id:
        mismatched_scan_mask.append(row['study_UID'])
    

if len(mismatched_scan_mask) > 0:
    for study_uid in mismatched_scan_mask:
        report.write(f"{study_uid} has a mismatched scan and mask. \n")
else: 
    report.write("No scans and masks are mismatched. \n")

report.close()