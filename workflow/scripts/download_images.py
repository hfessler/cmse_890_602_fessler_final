"""Download images from TCIA

    Download images from TCIA based on collection name, name of image, and patient ID.

    Args:
        collection: String of the name of the collection from TCIA. 
        patient_id: List of strings of patient IDs used in a collection.
        image_name: String detailing the type of image to be downloaded, as descried by the 
        collection.
        temp_storage: File path to a location to temporarly store downloads. 

    Saves each image to its own folder named via the patient ID and timing of the image. 

    Will rename folder containing download to be based on patient ID and time of image as 
    opposed to the default TCIA universal ID. While the universal IDs are unique, they are 
    very long and not human readable and contain no direct information about the image.   
    """

import os
from tcia_utils import nbia

for patient_id in snakemake.params.patient_ids:

    # Dataframe containing each image in the collection.
    df = nbia.getSeries(collection=snakemake.params.collection, format="df")

    # Filter this dataframe first by patient ID, then by image type.
    filtered_df = df.loc[df['PatientID'] == patient_id]
    filtered_df = filtered_df.loc[filtered_df['SeriesDescription']
                                  == snakemake.params.image_name]

    # uids are unique identifiers used by TCIA to identify every
    # image stored by them, and are required to download images.
    uids = filtered_df['SeriesInstanceUID'].tolist()
    nbia.downloadSeries(uids,
                        input_type="list",
                        path=snakemake.params.temp_storage)

    # uids are not connected to collection based ids
    # (like patient id or timing) so we will rename them to be more readable.
    for uid in uids:
        metadata = nbia.getSeriesMetadata(uid)[0]

        # ISPY2 images are uniquely names by the timing of the image
        # (found in study description) and the subject id.
        time = metadata['Study Description'][-1]
        subject_id = metadata['Subject ID']
        os.rename(os.path.join(snakemake.params.temp_storage, uid),
                  os.path.join(snakemake.params.temp_storage,
                               f'{subject_id}_{time}'))
