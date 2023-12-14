import pandas as pd
import os
from tcia_utils import nbia
from glob import glob

for patient_id in snakemake.params.patient_ids:
    df = nbia.getSeries(collection=snakemake.params.collection, format="df")
    filtered_df = df.loc[df['PatientID'] == patient_id]
    scans_df = filtered_df.loc[filtered_df['SeriesDescription'] == snakemake.params.scan]

    scan_uids = scans_df['SeriesInstanceUID'].tolist()

    nbia.downloadSeries(scan_uids, input_type = "list", path = snakemake.params.scan_storage)

    for uid in scan_uids:
        metadata = nbia.getSeriesMetadata(uid)[0]
        time = metadata['Study Description'][-1]
        subject_id = metadata['Subject ID']
        os.rename(os.path.join(snakemake.params.scan_storage,uid), os.path.join(snakemake.params.scan_storage,f'{subject_id}_{time}'))
