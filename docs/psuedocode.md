
The following is a high level psuedocode overview of the workflow. This walks through the whole process form downloading the images to moving the outputs to permanent storage.

```
IMPORT download_image FROM The Cancer Imaging Archieve Utilities

IMPORT convert_DICOM_to_NRRD, scan_processing, mask_processing FROM Simple Insight Toolkit

IMPORT extract_radiomics, get_bounding_box, crop_image FROM PyRadiomics

DEF download_images(patient_id, times_of_images, type_of_image)
    FOR time IN times_of_images:
        APPEND download_image(patient_id, time, type_of_image) TO images
    RETURN images

DEF process_images(images, type_of_image):
    FOR image IN images:
        IF type_of_image IS mask_type:
            APPEND mask_processing(image) TO processed_images
        IF type_of_image IS scan_type:
            APPEND scan_processing(image) TO processed_images
    FOR image IN processed_images:
        APPEND convert_DICOM_to_NRRD(image) TO NRRD_images
    RETURN NRRD_images

DEF create_output_files(patient_id, times_of_images)

    masks = download_images(patient_id, times_of_images, mask_type)
    masks_processed = process_images(masks, mask_type)

    scans = download_images(patient_id, times_of_images, scan_type)
    scans_processed = process_images(scans, scan_type)
    
    bounding_boxes = get_bounding_box(scans_processed, masks_processed)

    scans_cropped = crop_image(scans_processed, bounding_box)

    radiomics = extract_radiomics(scans_cropped, masks_processed)

    RETURN masks_processed, scans_processed, scans_cropped, bounding_boxes, radiomics

times = example_requested_times
FOR patient_id IN list_of_patient_ids:
    masks_processed, scans_processed, scans_cropped, bounding_boxes, radiomics = create_output_files(patient_id, times)

    MOVE masks_processed TO permanate_storage
    MOVE scans_processed TO permanate_storage
    MOVE scans_cropped TO permanate_storage
    MOVE bounding_boxes TO permanate_storage
    MOVE radiomics TO permanate_storage
```