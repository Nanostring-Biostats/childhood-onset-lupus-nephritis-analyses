This repo holds analyses used in the manuscript, "Single cell spatial transcriptomic profiling of childhood-onset lupus nephritis reveals complex interactions between kidney stroma and infiltrating immune cells".

Note: the datasets used here are too big to host on github. They can be downloaded from Figshare at: _____________________________________. 

How to use:

First, copy this repo to the machine where you plan to run the scripts. 
Next, download the datasets from the above url and place them in the copied repo's "processed_data" folder.
Finally, execute the R scripts in numerical order. They each execute in the directory in which they lie.

The cell typing script "3. cell typing by batch.R" uses the complete dataset "data_with_failed_slide.RData", which includes a failed slide omitted from downsteam analyses. 
All other analyses, and all results in the manuscript, use "cleaneddata.RData". 