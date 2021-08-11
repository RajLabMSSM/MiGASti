## Omics analysis for cultured and stimulated microglia 2nd pass

> This repository includes code and plots. Exploratory analysis and intermediate processing files are too large for this repository.

> This repository includes codes and plots after removal of some additional samples (RNA/DNA mismatch; sample swaps and donor 14-055 GFM has been changed to donor 14-051). See docs 1st_pass for first analyses, see docs 2nd_ pass for updated analyses.

1. [QC](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/QC_cor.html). Quality control of all samples including removal of additional samples.

2. Metadata
- [Metadata](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/Metadata_all.html). Organizing of the metadata and general demographics for the included samples after QC filtering of all samples. 
- [Metadata_cultured_stimulated](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/Metadata_cultured.html). Metadata and general demographics for the included samples after QC filtering of all cultured + stimulated samples.
- [Metadata_ununstim](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/Metadata_ununstim.html).
 
3.Exploratory plots: 
 - [Exploratory plots after filters](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210210_PCA_filtering.html). PCA's, heatmaps with linear regression between the first 20 PCs with  and without uncultured samples.  
 
4. Variance partition 
- [Variance partition_cultured](https://rajlabmssm.github.io/MiGASti/docs/2nd_passVariance_partition_cultured.html). Variance partition for only cultured samples with and without TNFa. Technical + biological factors combined and biological factors only. Technical factors was the collinearity to high, so could not be modeled.  

5. Internal QC; TNFa is already excluded. 
- [Expression of brain markers](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210604_Markers_braincells.html). Expression of microglia, astrocyte, oligodendrocyte, neuron markers in all samples.
- [Expression of myeloid markers](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210604_Markers_myeloid.html). Expression of monocyte, macrophage and microglia markers in all samples. 
- [Mitochondrial_genes](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210511_Mitochondrial_genes.html). Percentage of mitochondrial genes in total dataset after filtering out lowly expressed genes (< 1 in 50% of the samples).
- [Stimulations](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20211206_Stimulations_heatmap2.html). Heatmap of stimulation specific responses. Expression of ligands in vitro vs ex vivo.
- [Expression of AD genes](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210224_Markers_AD_genes.html). Boxplots showing expression of AD genes (based on TWAS AD and snRNAseq) in IL4 and DEX. 

6. DEG analysis (with only cultured samples)
- [DESeq2_GFM_all](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210223DiffExpression_GFM_all.html). DESeq2 analysis with contrasts of GFM samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_SVZ_all](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210223DiffExpression_SVZ_all.html). DESeq2 analysis with contrasts of SVZ samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_GTS_all](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210223DiffExpression_GTS_all.html). DESeq2 analysis with contrasts of GTS samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_CC_all](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210223DiffExpression_CC_all.html). DESeq2 analysis with contrasts of CC samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 
- [DESeq2_THA_all](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210225DiffExpression_THA_all.html). DESeq2 analysis with contrasts of THA samples only: stimuli vs unstim for all conditions. Number of differential expressed genes, Vulcano plots, MA plots, list of top genes. 

Note. SVZ showed most significant genes up/down after stimulation.

| GFM  	| FDR 5% 	| LogFC \|1\| 	| GTS 	| FDR 5% 	| logFC \|1\| 	| CC 	| FDR 5% 	| logFC \|1\| 	| THA 	| FDR 5% 	| logFC \|1\| 	| SVZ 	| FDR 5% 	| logFC \|1\| 	|
|------	|--------	|-------------	|-----	|--------	|-------------	|----	|--------	|-------------	|-----	|--------	|-------------	|-----	|--------	|-------------	|
| LPS  	| 376    	| 153         	|     	| 44     	| 18          	|    	| 90     	| 55          	|     	| 226    	| 95          	|     	| 1895   	| 450         	|
| IFNy 	| 155    	| 6           	|     	| 62     	| 59          	|    	| 76     	| 67          	|     	| 66     	| 61          	|     	| 292    	| 178         	|
| R848 	| 31     	| 13          	|     	| 0      	| 0           	|    	|        	|             	|     	|        	|             	|     	| 708    	| 176         	|
| DEX  	| 622    	| 542         	|     	|        	|             	|    	|        	|             	|     	|        	|             	|     	|        	|             	|
| IL4  	| 49     	| 47          	|     	|        	|             	|    	|        	|             	|     	|        	|             	|     	|        	|             	|
| ATP  	| 134    	| 134         	|     	|        	|             	|    	|        	|             	|     	|        	|             	|     	| 0      	| 0           	|

7. Differential expression across regions
- [DGE_SVZ_compared](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210303_DEG_FC_scatterplot.html). Scatterplots of logFC of genelist with all differential expressed genes compared between SVZ and the other region for LPS and IFNy only. 
- [DEG_SVZ](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210608_SVZXGFM_Limma_Deseq2.html). Comparison SVZxGFM in baseline samples (unstim) compared with DEG LPS SVZ + DEG IFNy SVZ + ex vivo SVZxGFM MiGA. DESeq2 used. 

8. DREAM only cultured + stimulated samples 
- [DREAM_analysis_cultured](https://rajlabmssm.github.io/MiGASti/docs/20212203_DREAM.html).
- [DREAM_analysis_uncultured](https://rajlabmssm.github.io/MiGASti/docs/20210609_DREAM_cultured_uncultured.html). 
- [DREAM_Vulcano](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210806_DREAM_volcanos_ms.html). Some vulcano plots for DEGs. 

12. Aging
- [Age_plots](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210806_age_plots_LPS_IFNy.html). Metadata, plots of effect of aging on immune response genes (6270 LPS and 79 IFNy genes) and overlapping genes between aging and immune response (151 LPS and 17 IFNy). Counts are voom normalized. Pairwise comparisons are not incorporated, since file is too big. 
- [Age_quantiles](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20211206_AGE_quantiles.html). MiGA aging genes as input to check expression of these genes after stimulation with LPS and IFNy or baseline. 

13. Differential transcript usage
- [DTU_heatmap](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210224_DTU_FC_heatmap_gene_names.html). Heatmap of significant DTU genes that overlap with DEG. 
- [DTU_plots_download](https://rajlabmssm.github.io/MiGASti/docs/20210511_DTU_plots_download.html). Summary of results including table for download of significant results. Code for DTU analysis can be found in docs 1st pass. Including boxplots showing direction of effects. 
- [DTU_overlaps_genes](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210512_DTU_genes_overlap.html). Overlap between transcripts DTU and genes per stimulation and significance of overlap between background of 20.000 genes. 

15. Comparison with monocyte data (Elisa) and microglia
- [Monocyte_comparison](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210806_DREAM_monocytes_comparison.html). Comparison of LPS and INFy genes up/down in microglia (6 hours) and monocytes (24 hours). microglia: logFC < -1 or > 1, monocytes: logFC < -5 or > 5. 

16. Enrichment analyses
Note. did not include TNFa, since stimulation did not work. Overlap between stimulations up/down FDR 5% and FDR 5% logFC < 0.5 or > 0.5. 

Input genelists:
Overlaps to highlight:
1. Genelist core microglia signature list n = 249 genes (Patir et al. 2019)

Functions:
2. GO: NFkB pathway n = 147 genes
3. GO: IFNy pathway n = 139 genes
4. GO: phagocytosis pathway n = X genes
5. Plaque induced genes n = 57 genes (Grubman et al. 2021)
6. Sensome genes n = 84 genes (Hickman et al. 2019)

Diseases:
7. snRNAseq Multiple sclerosis n = 64 genes (Schiermer et al. 2019), 
8. snRNAseq Alzheimers disesease n = 53 genes (Mathys et al. 2019), 
9. TWAS AD genes (Raj et al. 2018) 
10. TWAS PD genes (Li et al. 2019)
11. eQTL PD MiGA (Paiva de lopes et al)
12. eQTL AD MiGA (Paiva de Lopes et al)
13. eQTL MS MiGA (Paiva de Lopes et al)
15. New AD SNPs (Bellenguez et al.)
14. AD bulk brain. List of differential expressed genes in AD bulk brain compared to controls in 3 datasets (MAYO, ROSMAP, MSSM)
15. scRNAseq Human alzheimer genes up n = 22/down n = 53 genes + combined (Srinsivan et al. 2019) (HAM signature) 
16. scRNAseq of activation response microglia in AD mouse n = 108 genes (Sala-Frigerio et al. 2019) 6 scRNAseq disease associated microglia up n = 280 /down n = 22 genes (Keren-Shaul et al. 2019)

- [LPS_enrichment](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210624_genelist_enrichment_LPS_function_disease_FC0.5.Rmd)
- [IFNy_enrichment](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210624_genelist_enrichment_IFNy_function_disease_FC0.5.Rmd)
- [R848_enrichment](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210624_genelist_enrichment_R848_function_disease_FC0.5.Rmd)
- [IL4_enrichment](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210624_genelist_enrichment_IL4_function_disease_FC0.5.Rmd)
- [ATP_enrichment](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210624_genelist_enrichment_ATP_function_disease_FC0.5.Rmd)
- [DEX_enrichment](https://rajlabmssm.github.io/MiGASti/docs/2nd_pass/20210624_genelist_enrichment_DEX_function_disease_FC0.5.Rmd)

 
 
 
 
 




















