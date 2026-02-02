## UCSC Genome Browser Visualizations
Track data hubs can be used to visualize genomics signal tracks using the UCSC Genome Browser (for more information, please see the [Raney et al., 2014](https://academic.oup.com/bioinformatics/article/30/7/1003/232409) publication.

We developed a track data hub to visualize both maxATAC-derived _in-silico_ and experimental ChIP-seq datasets in tandem with chromatin accessibility data. Data can be explored in multiple ways using this tool. 

### Visualization of Genomic Regions
After loading the track data hub (link [here](https://gb.research.cchmc.org/hub/group/memoryCD4/hub.txt)) in "My Data" > "Track Hubs" on the [UCSC Genome Browser website](https://genome.ucsc.edu/index.html), you can select specific genomic locations to visualize all subpopulation-resolved maxATAC predictions and chromatin accessibility information in that area. For instance, the image below shows maxATAC predictions in each of the sixteen subpopulations at the promoter and gene body of _IFNG_.

<img width="2250" height="1612" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Screenshot" src="https://github.com/user-attachments/assets/f039fc1b-269c-4097-a52a-fa7e0f26d74b" />

### Visualization of Specific Subpopulations
If you would like to visualize maxATAC predictions for a subset of subpopulations (e.g., the Th1 resting and active subpopulations), navigate to the track settings page for this track and select the desired subopoulations to visualize, as shown below.

<img width="1482" height="882" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Settings" src="https://github.com/user-attachments/assets/62c322ae-6bb2-43dd-8147-2f20b08fe2c3" />

Once the changes are submitted, only the desired subpopulations appear in the Browser.

<img width="2250" height="441" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Th1_Screenshot" src="https://github.com/user-attachments/assets/7a7394f1-c465-4529-a4d0-675f63dbbce6" />

### Identification of Specific TFBS
By right-clicking on one of the subpopulation tracks and selecting the "pack" option, you can visualize all TFs with predicted binding sites in that region for the chosen subpopulation.

<img width="1225" height="605" alt="260201_CD4_T_Cell_UCSC_Genome_Browser_Settings_Two" src="https://github.com/user-attachments/assets/9184843d-9c1c-4db3-8c91-eaeffe032694" />

For example, the following screenshot demonstrates all of the TFs that are predicted to bind at one of the promoters of _IFNG_ in the resting Th1 subpopulation. 

![260201_CD4_T_Cell_UCSC_Genome_Browser_Th1_Rest_Screenshot](https://github.com/user-attachments/assets/386d0a5b-173b-4f49-9921-06ac818ab3be)


