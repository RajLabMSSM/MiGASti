## Omics analysis for cultured and stimulated microglia 

> This repository includes code and plots. Exploratory analysis and intermediate processing files are too large for this repository.

1. Metadata
- [Metadata](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Metadata.html). Organizing of the metadata and general demographics for the included samples after QC filtering of all samples. 
- [Metadata_cultured_stimulated](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Metadata_cultured.html). Metadata and general demographics for the included samples after QC filtering of all cultured + stimulated samples.
- [Metadata_ununstim](https://rajlabmssm.github.io/MiGASti/docs/tree/master/1st_pass/Metadata_ununstim.html).

2. [QC](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/QC_cor.html). Quality control of all samples.
 
3.Exploratory Analysis: 
 - [Exploratory plots before filters](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/20210209_PCA_heatmap_before_filtering.html). PCA's between the first 10 PCs and covariates with (533 samples) and without (483 samples) uncultured samples.
 - [Exploratory plots after filters](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210210_PCA_filtering.html). PCA's, heatmaps with linear regression between the first 20 PCs with (496 samples) and without uncultured (454 samples) samples.  
 
4. Variance partition 
- [Variance partition_all](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition.html). Variance partition for all samples.
- [Variance partition_cultured](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition_cultured.html). Variance partition for only cultured samples.
- [Variance partition_cultured_GFM](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition_GFM.html). Variance partition for only cultured samples.
- [Variance partition_cultured_GTS](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition_GTS.html). Variance partition for only cultured samples.
- [Variance partition_cultured_SVZ](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition_SVZ.html). Variance partition for only cultured samples.
- [Variance partition_cultured_THA](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition_THA.html). Variance partition for only cultured samples.
- [Variance partition_cultured_CC](https://rajlabmssm.github.io/MiGASti/tree/master/docs/1st_pass/Variance_partition_CC.html). Variance partition for only cultured samples.

5. Internal QC
- [Expression of pre-selected markers specific for culturing](https://rajlabmssm.github.io/MiGASti/tree/master/docs/20210217_Markers_homeostatic.html). Boxplots with TPM expression of genes that are expected to go down after culturing based on the literature.
- [Expression of pre-selected markers specific for stimulated conditions](https://rajlabmssm.github.io/MiGASti/tree/masterdocs/20210217_Markers_allstims.html). Boxplots with TPM expression for genes that respond to specific stimuli and heatmap with all markers combined. 
- [Expression of pre-selected markers specific for TNFa stimulation](https://rajlabmssm.github.io/MiGASti/tree/master/docs/20210225_Markers_TNFa.html). Boxplots with genes that respond to TNFa in cultured monocytes (effect most striking after 24h).
- [Expression of pre-selected markers for apoptosis/cell death](https://rajlabmssm.github.io/MiGASti/tree/master/docs/120210224_Markers_apoptotic.html). Boxplots with TPM expression of genes that are involved in apoptotic processes (CASP3; specific) 
- [Expression of pre-selected neurotransmitter markers](https://rajlabmssm.github.io/MiGASti/docs/20210304_Markers_neurotransmitters.html). Boxplots with TPM expression of genes that are involved neurotransmitter activity. 
- [Expression of neurodevelopmental markers](https://rajlabmssm.github.io/MiGASti/docs/20211206_Lot_genes.html). 
- [Expression of brain markers](https://rajlabmssm.github.io/MiGASti/docs/20210604_Markers_braincells.html). Expression of microglia, astrocyte, oligodendrocyte, neuron markers in all samples
- [Expression of myeloid markers](https://rajlabmssm.github.io/MiGASti/docs/20210604_Markers_myeloid.html). Expression of monocyte, macrophage and microglia markers in all samples. 
- [Mitochondrial_genes](https://rajlabmssm.github.io/MiGASti/docs/20210511_Mitochondrial_genes.html). Percentage of mitochondrial genes in total dataset after filtering out lowly expressed genes (< 1 in 50% of the samples).
- [Signaling_pathways](https://rajlabmssm.github.io/MiGASti/docs/20210604_Signaling_pathways.html). Heatmap of expression of ligands for different stimulations ex vivo vs in vitro. 
- [Stimulations](https://rajlabmssm.github.io/MiGASti/docs/20211206_Stimulations_heatmap2.html). Heatmap of stimulation specific responses. Expression of ligands in vitro vs ex vivo.
- [Mismatch_RNAseq](https://rajlabmssm.github.io/MiGASti/docs/20210606_sex_mismatch_RNA.html). Sample swaps RNA seq based on reported sex. 
- [Samples_brain_markers](https://rajlabmssm.github.io/MiGASti/docs/20210604_Markers_brain_sample_check.html). Individual check of samples on brain markers. Exclusion of MG-01-RNA, MG-02-RNA and MG-03-RNA.
- [Samples_myeloid_markers](https://rajlabmssm.github.io/MiGASti/docs/20210604_Markers_brain_sample_check_myeloid.html). Individual check of samples on myeloid markers. Exclusion of MG-01-RNA, MG-02-RNA and MG-03-RNA.

6. DEG analysis (with only cultured samples)
- [DESeq2_GFM_subset](https://rajlabmssm.github.io/MiGASti/docs/20210217_DiffExpression_GFM.html). DESeq2 analysis of GFM samples: unstimulated samples (baseline) compared to LPS and IFNy stimulation seperatly. 
- [DESeq2_GFM_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_GFM_all.html). DESeq2 analysis with contrasts of GFM samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_SVZ_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_SVZ_all.html). DESeq2 analysis with contrasts of SVZ samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_GTS_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_GTS_all.html). DESeq2 analysis with contrasts of GTS samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_CC_all](https://rajlabmssm.github.io/MiGASti/docs/20210223DiffExpression_CC_all.html). DESeq2 analysis with contrasts of CC samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_THA_all](https://rajlabmssm.github.io/MiGASti/docs/20210225DiffExpression_THA_all.html). DESeq2 analysis with contrasts of THA samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 

7. Differential expression across regions
- [DGE_all_stimuli_across_regions](https://rajlabmssm.github.io/MiGASti/docs/20210224_DEG_FC_heatmap_gene_names.html). Heatmaps, PCAs, upset plots with differentially expressed genes per stimulation seperate across brain regions with filter Log FC > 1 or Log FC < -1. 
- [DGE_2_regions_compared](https://rajlabmssm.github.io/MiGASti/docs/20210303_DEG_FC_scatterplot.html). Scatterplots of logFC of genelist with all differential expressed genes compared between two brain regions for a subset of the different stimuli (LPS, IFNy, R848, TNFa). No logFC treshold.

8. DREAM only cultured + stimulated samples 
- [DREAM_analysis](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Contrast plot. 
- [DREAM_plots](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_plots_latest.html). Number of differential expressed genes, Vulcano plots, MA plots, FDR distribution plots, FC plots. 
- [DREAM_boxplots](https://rajlabmssm.github.io/MiGASti/docs/20210331_DREAM_Boxplots_tpm.html). Directionality of top 6 genes for each condition. Genes counts + log2((tpm)+1).
- [DREAM_lists](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM_DEG_download.html).  DEG genelists for download. No logFC treshold. 
- [DREAM_downstream_all](https://rajlabmssm.github.io/MiGASti/docs/20210514_DEG_FC_heatmap_DREAM_tresholds.html). Heatmaps, PCAs, upsetplot, scatterplots showing results of different stimuli for FDR 5% log FC > 1 or < -1 (1204 genes), FDR 5% in 2 stimulations (3161 genes) and FDR 5% (8521 genes). 
- [DREAM_downstream_LPS_IFNy](https://rajlabmssm.github.io/MiGASti/docs/20210514_DEG_FC_heatmap_DREAM_tresholds_LPS_IFNy.html), Heatmaps, PCAs, upsetplot, scatterplots showing results of LPS and IFNy for FDR 5% log FC > 1 or < -1 (210 genes), FDR 5% in 2 stimulations (348) and FDR 5% (6571).


9. K-means clustering on all cultured + stimulated samples
- [kmeans_clustering_plots_all](https://rajlabmssm.github.io/MiGASti/docs/20210514_kmeans_dream_stimulations.html). Elbow, silhouette and gap statistic method for kmeans clustering of genes and some trail and error with clustering settings with use of gprofiler for FDR 5% log FC > 1 or < -1 (1204 genes).

10. K-means clustering on LPS and IFNy samples
- [kmeans_clustering_plots_LPS_IFNy_FDR5](https://rajlabmssm.github.io/MiGASti/docs/2021049_kmeans_dream_LPS_IFNy_FDR5.html). Elbow, silhouette and gap statistic method for kmeans clustering of genes and some trail and error with clustering settings with use of gprofiler for FDR 5% (6571) of a subset of LPS and IFNy genes. 
- [kmeans_clustering_plots_LPS_IFNy_LOGFC](https://rajlabmssm.github.io/MiGASti/docs/2021049_kmeans_dream_LPS_IFNy_logFC1.html). Elbow, silhouette and gap statistic method for kmeans clustering of genes and some trail and error with clustering settings with use of gprofiler for logFC > 1 or logFC < -1, FDR 5% (210) of a subset of LPS and IFNy genes.
- [kmeans_heatmap](https://rajlabmssm.github.io/MiGASti/docs/2021049_heatmap_kmeans.html). Heatmap kmeans clustering with 4 clusters.

11. DREAM all samples
- [DREAM_analysis](https://rajlabmssm.github.io/MiGASti/docs/20210609_DREAM_cultured_uncultured.html). Cultured versus uncultured. Not corrected for region, since 3 regions were NA. Dream model. Number of differential expressed genes, Volcano plots, MA plots, FDR distribution plots, FC plots, boxplots top 6 genes and table for download. 
- [DREAM_Vulcano](https://rajlabmssm.github.io/MiGASti/docs/20210806_DREAM_volcanos_ms.html). Some vulcano plots for DEGs. 

12. Aging
- [DREAM_age](https://rajlabmssm.github.io/MiGASti/docs/20210806_DREAM_AGE.html). Age as coefficient in all cultured and stimulated samples. Number of differential expressed genes, Volcano plots, MA plots, FDR distribution plots, FC plots.  
- [Age_plots](https://rajlabmssm.github.io/MiGASti/docs/20210806_age_plots_LPS_IFNy.html). Metadata, plots of effect of aging on immune response genes (6270 LPS and 79 IFNy genes) and overlapping genes between aging and immune response (151 LPS and 17 IFNy). Counts are voom normalized. Pairwise comparisons are not incorporated, since file is too big. 
- [DREAM_age_stimulation](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html). Interaction analysis age * stimulation. 0 differentially expressed genes in LPS and IFNy associated with aging. 

13. Differential transcript usage
- [DTU_plots_download](https://rajlabmssm.github.io/MiGASti/docs/20210511_DTU_plots_download.html). Summary of results including table for download of significant results. Code for DTU analysis can be found in docs. 
- [DTU_overlaps_genes](https://rajlabmssm.github.io/MiGASti/docs/20210512_DTU_genes_overlap.html). Overlap between transcripts DTU and genes per stimulation and significance of overlap between background of 20.000 genes. 
- [DTU_Vulcano](https://rajlabmssm.github.io/MiGASti/docs/20210806_DTU_volcanos_ms.html). Some vulcano plots for DTU. 

14. Sex
- [DREAM_sex](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20211305_DREAM_SEX_CHRXY.html). Differential expression analysis male/female in all cultured samples. 0 differentially expressed genes.

15. Comparison with monocyte data (Elisa) and microglia (MiGA), macrophages (Yang) and monocytes/microglia (Gosselin) 
- [Ranking_myeloid_genes](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210514_monocyte_microglia_ranking_top400.html). Ranking of microglia/monocyte/macrophage genes across the myeloid datasets top 400. 
- [PCA_myeloid](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210514_monocyte_microglia_macrophage.html). PCA plot of different myeloid datasets (MiGA, Gosselin, monocytes, macrophages).
- [PCA_stimulation](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210514_monocyte_microglia.html). PCA plot of cultured LPS and IFNy microglia, cultured, LPS and IFNy monocytes and baseline microglia.
- [Monocyte_comparison](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210806_DREAM_monocytes_comparison.html). Comparison of LPS and INFy genes up/down in microglia (6 hours) and monocytes (24 hours). microglia: logFC < -1 or > 1, monocytes: logFC < -5 or > 5. 

16. Enrichment analyses
Note. did not include TNFa, since stimulation did not work. Overlap between stimulations up/down FDR 5% and FDR 5% logFC < -1 or > 1. 

Input genelists:
1. Genelist core microglia signature list n = 249 genes (Patir et al. 2019)
2. snRNAseq Multiple sclerosis n = 64 genes (Schiermer et al. 2019), 
3. snRNAseq Alzheimers disesease n = 53 genes (Mathys et al. 2019), 
4. mouse snRNAseq response to demyelinatin injury after lysolecithin injection into the brain n = 119 genes (Hammond et al. 2019)  
5. scRNAseq of activation response microglia in AD mouse n = 108 genes (Sala-Frigerio et al. 2019)
6 scRNAseq disease associated microglia up n = 280 /down n = 22 genes (Keren-Shaul et al. 2019), 
7. scRNAseq microglia from mouse stimulated with LPS 3 hours n = 367 genes (acute LPS; Gerrits et al. 2020),
8. scRNAseq microglia from mouse stimulated with LPS 24 hours n = 152 genes (subacute LPS; Sousa et al. 2020), 
9. GO: NFkB pathway n = 147 genes
10. GO: IFNy pathway n = 139 genes
11. GO: cold stress response pathway n = 42 genes
12. Plaque induced genes n = 57 genes (Grubman et al. 2021)
13. Sensome genes n = 84 genes (Hickman et al. 2019)
14. TWAS AD genes (Raj et al. 2018) 
15. TWAS PD genes (Li et al. 2019)
16. scRNAseq Human alzheimer genes up n = 22/down n = 53 genes (Srinsivan et al. 2019)
17. Homeostatic genes up/down n = 44 genes in the mature homeostatic microglia phenotype based on bulk RNAseq Ormel et al. 2018, Gosselin et al. 2017 and Kracht et al. 2020).
18. bulk RNAseq cultured microglia compared to ex vivo up n = 3674 / down n = 4121 genes (Gosselin et al. 2017)
19. Immediate early response genes (IEGs), genes whose transcripts are rapidly expressed after cellular stimulation, without the need for new protein synthesis to control this expression n = 128 genes (Wu et al. 2017) 
20. Microglia genes up n = 338/down n = 1355 in aging (MIGA)

- [LPS_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_LPS.html)
- [IFNy_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_IFNy.html)
- [R848_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_R848.html)
- [IL4_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_IL4.html)
- [ATP_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_ATP.html)
- [DEX_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_DEX.html)
- [Cluster_LPS_IFNy_enrichment](https://rajlabmssm.github.io/MiGASti/blob/master/docs/1st_pass/20210624_genelist_enrichment_clusters.html). Geneset enrichment for clusters. 