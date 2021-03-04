## Omics analysis for cultured and stimulated microglia 

> This repository includes code and plots. Exploratory analysis and intermediate processing files are too large for this repository.

1. [Metadata](https://rajlabmssm.github.io/MiGASti/docs/Metadata.html). Organizing of the metadata and general demographics for the included samples after QC filtering.

2. [QC](https://rajlabmssm.github.io/MiGASti/docs/QC_cor.html). Quality control of all samples.
 
3.Exploratory Analysis: 
 - [Exploratory plots before filters](https://rajlabmssm.github.io/MiGASti/docs/20210209_PCA_heatmap_before_filtering.html). PCA's between the first 10 PCs and covariates with (533 samples) and without (483 samples) uncultured samples.
 - [Exploratory plots after filters](https://rajlabmssm.github.io/MiGASti/docs/20210210_PCA_filtering.html). PCA's, heatmaps with linear regression between the first 20 PCs with (496 samples) and without uncultured (454 samples) samples.  
 
4. [Variance partition](https://rajlabmssm.github.io/MiGASti/docs/Variance_partition.html) 

5. Expression of stimulation specific markers
- [Expression of pre-selected markers specific for culturing](https://rajlabmssm.github.io/MiGASti/docs/20210217_Markers_homeostatic.html). Boxplots with TPM expression of genes that are expected to go down after culturing based on the literature.
- [Expression of pre-selected markers specific for stimulated conditions](https://rajlabmssm.github.io/MiGASti/docs/20210217_Markers_allstims.html). Boxplots with TPM expression for genes that respond to specific stimuli and heatmap with all markers combined. 
- [Expression of pre-selected markers specific for TNFa stimulation](https://rajlabmssm.github.io/MiGASti/docs/20210225_Markers_TNFa.html). Boxplots with genes that respond to TNFa in cultured monocytes (effect most striking after 24h).
- [Expression of pre-selected markers for apoptosis/cell death](https://rajlabmssm.github.io/MiGASti/docs/20210224_Markers_apoptotic.html). Boxplots with TPM expression of genes that are involved in apoptotic processes. 
- [Expression of pre-selected neurotransmitter markers](https://rajlabmssm.github.io/MiGASti/docs/20210304_Markers_neurotransmitters.html). Boxplots with TPM expression of genes that are involved neurotransmitter activity. 

6. DEG analysis
- [DESeq2_GFM_subset](https://rajlabmssm.github.io/MiGASti/docs/20210217_DiffExpression_GFM.html). DESeq2 analysis of GFM samples: unstimulated samples (baseline) compared to LPS and IFNy stimulation seperatly. 
- [DESeq2_GFM_all]((https://rajlabmssm.github.io/MiGASti/docs/20210223_DiffExpression_GFM_all.html). DESeq2 analysis with contrasts of GFM samples only: stimuli vs unstim for all conditions. Number of differentially expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_SVZ_all]((https://rajlabmssm.github.io/MiGASti/docs/20210223_DiffExpression_SVZ_all.html). DESeq2 analysis with contrasts of SVZ samples only: stimuli vs unstim for all conditions. Number of differentially expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_GTS_all]((https://rajlabmssm.github.io/MiGASti/docs/20210223_DiffExpression_GTS_all.html). DESeq2 analysis with contrasts of GTS samples only: stimuli vs unstim for all conditions. Number of differentially expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_CC_all]((https://rajlabmssm.github.io/MiGASti/docs/20210223_DiffExpression_CC_all.html). DESeq2 analysis with contrasts of CC samples only: stimuli vs unstim for all conditions. Number of differentially expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_THA_all]((https://rajlabmssm.github.io/MiGASti/docs/20210223_DiffExpression_THA_all.html). DESeq2 analysis with contrasts of THA samples only: stimuli vs unstim for all conditions. Number of differentially expressed genes, Vulcano plots, MA plots, list of top genes. 

7. Differential expression across regions
- [DGE_all_stimuli_across_regions]((https://rajlabmssm.github.io/MiGASti/docs/20210224_DEG_FC_heatmap.html). Heatmaps, PCAs, upset plots with differentially expressed genes per stimulation seperate across brain regions with filter Log FC > 1 or Log FC < -1. 
- [DGE_2_regions_compared](https://rajlabmssm.github.io/MiGASti/docs/20210223_DEG_FC_scatterplot.html). Scatterplots of logFC of genelist with all differentially expressed genes compared between two brain regions for a subset of the different stimuli (LPS, IFNy, R848, TNFa)

8. DREAM





