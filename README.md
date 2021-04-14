## Omics analysis for cultured and stimulated microglia 

> This repository includes code and plots. Exploratory analysis and intermediate processing files are too large for this repository.

1. Metadata
- [Metadata](https://rajlabmssm.github.io/MiGASti/docs/Metadata.html). Organizing of the metadata and general demographics for the included samples after QC filtering of all samples. 
- [Metadata_cultured_stimulated](https://rajlabmssm.github.io/MiGASti/docs/Metadata_cultured.html). Metadata and general demographics for the included samples after QC filtering of all cultured + stimulated samples.
- [Metadata_ununstim](https://rajlabmssm.github.io/MiGASti/docs/Metadata_ununstim.html).

2. [QC](https://rajlabmssm.github.io/MiGASti/docs/QC_cor.html). Quality control of all samples.
 
3.Exploratory Analysis: 
 - [Exploratory plots before filters](https://rajlabmssm.github.io/MiGASti/docs/20210209_PCA_heatmap_before_filtering.html). PCA's between the first 10 PCs and covariates with (533 samples) and without (483 samples) uncultured samples.
 - [Exploratory plots after filters](https://rajlabmssm.github.io/MiGASti/docs/20210210_PCA_filtering.html). PCA's, heatmaps with linear regression between the first 20 PCs with (496 samples) and without uncultured (454 samples) samples.  
 
4. Variance partition 
- [Variance partition_all](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for all samples.
- [Variance partition_cultured](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for only cultured samples.
- [Variance partition_cultured_GFM](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for only cultured samples.
- [Variance partition_cultured_GTS](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for only cultured samples.
- [Variance partition_cultured_SVZ](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for only cultured samples.
- [Variance partition_cultured_THA](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for only cultured samples.
- [Variance partition_cultured_CC](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html). Variance partition for only cultured samples.

5. Expression of stimulation specific markers
- [Expression of pre-selected markers specific for culturing](https://rajlabmssm.github.io/MiGASti/docs/20210217_Markers_homeostatic.html). Boxplots with TPM expression of genes that are expected to go down after culturing based on the literature.
- [Expression of pre-selected markers specific for stimulated conditions](https://rajlabmssm.github.io/MiGASti/docs/20210217_Markers_allstims.html). Boxplots with TPM expression for genes that respond to specific stimuli and heatmap with all markers combined. 
- [Expression of pre-selected markers specific for TNFa stimulation](https://rajlabmssm.github.io/MiGASti/docs/20210225_Markers_TNFa.html). Boxplots with genes that respond to TNFa in cultured monocytes (effect most striking after 24h).
- [Expression of pre-selected markers for apoptosis/cell death](https://rajlabmssm.github.io/MiGASti/docs/20210224_Markers_apoptotic.html). Boxplots with TPM expression of genes that are involved in apoptotic processes. 
- [Expression of pre-selected neurotransmitter markers](https://rajlabmssm.github.io/MiGASti/docs/20210304_Markers_neurotransmitters.html). Boxplots with TPM expression of genes that are involved neurotransmitter activity. 
- [Expression of microglia markers](https://rajlabmssm.github.io/MiGASti/docs/20210604_Markers_braincells.html). Expression of microglia, astrocyte, oligodendrocyte, neuron markers in all samples
- [Expression of myeloid markers](https://rajlabmssm.github.io/MiGASti/docs/20210604_Markers_myeloid.html). Expression of monocyte, macrophage and microglia markers in all samples. 


6. DEG analysis (with only cultured samples)
- [DESeq2_GFM_subset](https://rajlabmssm.github.io/MiGASti/docs/20210217_DiffExpression_GFM.html). DESeq2 analysis of GFM samples: unstimulated samples (baseline) compared to LPS and IFNy stimulation seperatly. 
- [DESeq2_GFM_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_GFM_all.html). DESeq2 analysis with contrasts of GFM samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_SVZ_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_SVZ_all.html). DESeq2 analysis with contrasts of SVZ samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_GTS_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_GTS_all.html). DESeq2 analysis with contrasts of GTS samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_CC_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_CC_all.html). DESeq2 analysis with contrasts of CC samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_THA_all](https://rajlabmssm.github.io/MiGASti/docs/20210225DiffExpression_THA_all.html). DESeq2 analysis with contrasts of THA samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 

7. Differential expression across regions
- [DGE_all_stimuli_across_regions](https://rajlabmssm.github.io/MiGASti/docs/20210224_DEG_FC_heatmap_genes_names.html). Heatmaps, PCAs, upset plots with differentially expressed genes per stimulation seperate across brain regions with filter Log FC > 1 or Log FC < -1. 
- [DGE_2_regions_compared](https://rajlabmssm.github.io/MiGASti/docs/20210303_DEG_FC_scatterplot.html). Scatterplots of logFC of genelist with all differential expressed genes compared between two brain regions for a subset of the different stimuli (LPS, IFNy, R848, TNFa). No logFC treshold.

8. DREAM only cultured + stimulated samples 
- [DREAM_analysis](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Contrast plot. 
- [DREAM_plots](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Number of differential expressed genes, Vulcano plots, MA plots, FDR distribution plots, FC plots. 
- [DREAM_boxplots](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_Boxplots.html). Directionality of top 6 genes for each condition.
- [DREAM_lists](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_DEG_download.html).  DEG genelists for download. No logFC treshold. 
- [DREAM_downstream_logFC](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Heatmaps, PCAs, upsetplot, scatterplots showing results of different stimuli for logFC > 1 or -1. 
- [DREAM_downstream_FDR5](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Heatmaps, PCAs, upsetplot, showing results of different stimuli for FDR < 0.05.
- [DREAM_downstream_FDR5_2stim](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Heatmaps, PCAs, upsetplot, showing results of different stimuli for FDR < 0.05 in 2 or more stimuli.

9. K-means clustering on cultured + stimulated samples
- [kmeans_clustering_plots_logFC](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Elbow, silhouette and gap statistic method for kmeans clustering of genes logFC > 1 or -1. 
- [kmeans_clustering_plots_FDR5](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Elbow, silhouette and gap statistic method for kmeans clustering of genes FDR < 5.
- [kmeans_clustering_plots_FDR5_2stim](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots.html). Elbow, silhouette and gap statistic method for kmeans clustering of genes FDR < 5 in 2 or more stimulations.
 
10. DREAM all samples
- [DREAM_analysis](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Contrast plot.
- [DREAM_analysis](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). DEG table for download. 

11. Aging 
- [DREAM_age](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Age as coefficient in all cultured and stimulated samples. 
- [DREAM_age_plots](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Number of differential expressed genes, Vulcano plots, MA plots, FDR distribution plots, FC plots.  
- [DREAM_age_stimulation](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Interaction analysis age * stimulation. Differentially expressed genes per stimulation associated with aging. 

12. Differential transcript usage
- [DTU](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). 