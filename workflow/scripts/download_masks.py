import pandas as pd
import os
from tcia_utils import nbia
from glob import glob

for patient_id in snakemake.params.patient_ids:
    df = nbia.getSeries(collection=snakemake.params.collection, format="df")
    filtered_df = df.loc[df['PatientID'] == patient_id]
    masks_df = filtered_df.loc[filtered_df['SeriesDescription']
                               == snakemake.params.mask_name]

    mask_uids = masks_df['SeriesInstanceUID'].tolist()

    nbia.downloadSeries(mask_uids, input_type="list",
                        path=snakemake.params.mask_storage)

    for uid in mask_uids:
        metadata = nbia.getSeriesMetadata(uid)[0]
        time = metadata['Study Description'][-1]
        subject_id = metadata['Subject ID']
        os.rename(os.path.join(snakemake.params.mask_storage, uid), os.path.join(
            snakemake.params.mask_storage, f'{subject_id}_{time}'))
