import pandas as pd
from tcia_utils import nbia

df = nbia.getSeries(collection=snakemake.params.collection, format="df")
print(df['PatientID'])
filtered_df = df.loc[df['PatientID'] == snakemake.params.patient_id]
scans_df = filtered_df.loc[filtered_df['SeriesDescription'] == snakemake.params.scan]
masks_df = filtered_df.loc[filtered_df['SeriesDescription'] == snakemake.params.mask]

scan_uids = scans_df['SeriesInstanceUID'].tolist()
mask_uids = scans_df['SeriesInstanceUID'].tolist()

nbia.downloadSeries(scan_uids, input_type = "list", path = snakemake.params.scan_storage)
nbia.downloadSeries(mask_uids, input_type = "list", path = snakemake.params.mask_storage)
