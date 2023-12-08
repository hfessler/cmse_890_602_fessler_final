import pandas as pd
from tcia_utils import nbia

df = nbia.getSeries(collection=snakemake.params.collection, format="df")

filtered_df = df.loc[df['SeriesDescription'] == snakemake.params.scan]

output_df = pd.DataFrame(columns=['study_ID', 'scan_UID', 'mask_UID'])
output_df.to_csv(snakemake.output.csv)

scan_UIDs = df['SeriesInstanceUID'].tolist()
output_data = []
for scan_UID in scan_UIDs:
    # Each study (imaging event) is assigned a unique ID
    study_UID = df.loc[df['SeriesInstanceUID']
                       == scan_UID]['StudyInstanceUID'].item()
    # Limit the dataframe to only those of the study for the current scan
    study_df = df.loc[df['StudyInstanceUID'] == study_UID]

    # Find the mask UID associated with this study 
    mask_df = study_df[study_df['SeriesDescription']
                        == snakemake.params.mask]['SeriesInstanceUID']

    # Not every study produces a mask. If so, mark no mask found. 
    if mask_df.empty:
        mask_UID = pd.NA
    else:
        mask_UID = mask_df.item()
    output_data.append({'study_UID': study_UID, 'scan_UID': scan_UID, 'mask_UID': mask_UID})

output_df = pd.DataFrame(output_data)

output_df.to_csv(snakemake.output.csv, index=False)
