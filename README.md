# AMPSEQ

Combinatorial indexing during library preparation tags PCR products in an Illumina Plate with unique, per-well inline barcodes pairs. By matching reads to their expected inline barcodes, indigenous PCR products in a well can be differentiated from alien reads that intruded in the well during PCR setup.

The Amplicon Decontamination workflow in AmpSeq leverages this principle for detecting, reporting, and filtering out contamination. When using inline barcodes, AmpSeq seamlessly pipes decontamined read files into the Denoising workflow (if no barcodes are used, the workflow will assume reads have been correctly imputed to their well and automatically proceed with the next step, denoising).

Amplicon Denoising is an algorithmic process for distinguishing true base calls in amplicon sequencing data. Denoising methods identify and remove PCR errors introduced during amplification, errors in sequence determination during sequencing reactions, and erroneuos sequences, such as chimeras. It is a necessary step to enhance the fidelity of reported haplotypes, allowing for reliable downstream analysis and interpretations. To perform the denoising of reads, AmpSeq on Terra incorporates the DADA2 algorithm.

# USAGE 

## Uploading data

AmpSeq has been optimized to run in [Terra](https://app.terra.bio/). Running the pipeline requires creating a New Folder in the "DATA" tab of the user's Terra workspace. For this, the user must click the "DATA" tab and then the "Files" tab on the bottom left. After accessing the files section, the user must click the "New Folder" button on the top right. This will redirect the user to the new folder. Inside the new folder, the user can upload the data described in the previous section (The workflow will dynamically identify the files using their names. This is the reason why naming conventions must be followed strictly). After the data have been uploaded, the folder must be "staged" by importing the data as described in the Terra documentation.

## Running the workflow

Once inside the workflow, the user must select the "plate" in which he wants to perform the analyses. This is done through the "Select Data" button. In the path_to_fq field, the user must type "this.study_id". The remaining parameters are preconfigured and sufficient to fulfill the requirements of typical analyses. Nonetheless, users retain the flexibility to adjust these parameters to suit specific objectives. The extended documentation - link pending provides information about this configurable parameters. Finally, the user must launch the workflow by clicking the "SAVE" and "RUN ANALYSIS" buttons.

## Contact

Jorge Eduardo Amaya Romero:

jamayaro@broadinstitute.org
