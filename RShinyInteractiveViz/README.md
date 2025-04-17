## Rshiny Interactive Visualization
The RShiny application allows interactive visualization of gene expression, transcription factor activity, data dimension reduction, and more!
#### Line plot
The first tab in the RShiny app allows for visualization of gene expression over the T-cell activation timecourse. Begin by loading in the pre-processed RNA counts matrix and meta file stored inside the RShiny folder. Typing in any gene of interest will automatically display the lineplot for these gene for all of the celltypes in the dataset
<img width="600" alt="Screen Shot 2025-03-28 at 9 56 33 AM" src="https://github.com/user-attachments/assets/b734ea25-7402-4f1f-b926-913e8a5ab153" />
#### UMAP plot
The second tab is designed for visualization of the single-cell dimension reductions. First, load in the seurat object that has been optimized to load quickly into the app. Next, you can select from the 3 available types of dimension reduction available: Weighted Nearest Neighbor (WNN), RNA UMAP, and ATAC UMAP. All three of these reductions have undergone harmony integration to remove batch effect from the visualization. You can then type any gene of interest to view the expression pattern on the dimension reduction
<img width="600" alt="Screen Shot 2025-03-28 at 9 59 04 AM" src="https://github.com/user-attachments/assets/350cb5c4-29da-41c8-8483-2b80240a3ba1" />
#### Heatmap
The third tab allows for customizable heatmaps for visualization of both gene expression and estimated transcription factor activity (TFA) from our gene regulatory networks. Begin by loading in either the pre-processed RNA counts or activity matrix. Loading in these matrices should not take more than 30 seconds of processing. 
You must also load in the metadata file, which will be the same regardless if you are visualizing expression or activity. You can then select the celltypes you want to include in the heatmap using the checkboxes. Finally, enter a comma separated list of genes you want to visualize. The heatmaps have been programmed to automatically resize based on the number of genes requested, ensuring the size of each row is always optimized. 
<img width="600" alt="Screen Shot 2025-03-28 at 10 03 12 AM" src="https://github.com/user-attachments/assets/ddf85586-e9eb-4e3f-9372-1c618bf6ccc1" />
<img width="600" alt="Screen Shot 2025-03-28 at 10 06 00 AM" src="https://github.com/user-attachments/assets/3692aaa3-f503-493e-9625-20177b481b1e" />
