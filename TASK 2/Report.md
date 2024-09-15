<!--StartFragment-->

 <h1>PROGRAMMING IN BIOLOGY & DATA VISUALIZATION</h1>

   **“Analysis of gene expression raw data to predict differentially expressed genes of significance and perform functional enrichment analysis”**

**TEAM:**\
Ghizal Niko (@Ghizal), David Ojo (@Dajom), Ademilola Adekoya (@Demi\_Nevaeh), Eman Binte Hafeez (@emanuts), Lillian Mwinja (@LillianMwinja), Irina Andriushchenko (@rirkella), Antara Ghanta (@Tara) and Narendran Mayilvaganan (@Naren037)

**CODE:**  https://github.com/Naren037/hackbio-cancer-internship/blob/91f8a3e62d47e7301b609b618edb5a9efb4fe894/TASK%202/Code.R    

1. **INTRODUCTION**

   Heatmaps are graphical data representations using colours to indicate different values, with hierarchical clustering grouping rows or columns by similarity to facilitate pattern identification and analysis (Gu et al., 2016). In glioblastoma       research, heatmaps help identify molecular subtypes and gene signatures (Clough & Barrett, 2016) aiding biomarker discovery and targeted therapy (Verhaak et al., 2010). Fold change measures the difference in gene or protein expression levels between two conditions, and the p-value (usually < 0.05), which assesses the significance of these differences in hypothesis testing, is crucial in gene expression analysis. Functional enrichment analysis further identifies whether the genes with significant expression changes are overrepresented in specific biological pathways or functions. Combining these techniques allows researchers to effectively navigate the complexities of gene expression data and advance our understanding of biological processes.

2. **METHODOLOGY**

   528 differentially expressed genes from 10 glioblastoma samples were analyzed using RStudio. Heatmaps were generated using heatmap.2 and were visualised with sequential and diverging colour palettes. Clustering was performed by rows, columns, and both. Fold change was calculated using log2(mean\_group1) - log2(mean\_group2), and the p-value was calculated using paired t-tests. The upregulated and downregulated genes were subset using a fold change threshold of 1.5 and p-value < 0.05. The biological pathways involved were identified using ShinyGo 8.0 with an FDR cut-off of 0.16. A bubble plot was generated for the top 6 pathways. 

3. **RESULTS**

   3.1. **Heatmap analysis:**

      Two distinct groups were identified by visualising the patterns when clustered by columns. Samples S1 to S5 were taken as group 1, and S6 to S10 were as group 2.

      ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXekBQMwssEnQR7MvzeiRu0jTKhi9h5Rd5ZoeZz61NtqxUviTsFA5qXImCMG9TqNNBdKjwvTcSipwV4ew8A2orjLQLa4jZfOx3Dhl6S4Ul8ySZKFlRrmqzPkSImtVTVotaZRklIIaLWnzOnzXIhGZ5bANOY8?key=1GDaguhNYDgftMhdZ7Y7yw)

   **Figure 1.** Heat map generated using a sequential color palette and clustered by columns

      ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXduBfqOXAyKinYh4cKrzOD00gXDISbizNiyW8eDbiLW0RF9pN6hVnFzq6ZB6zjL149DFJ2ce5Kc9R5N6NSETXEoF1ETDWajiBEysfmxPr-yfU7721a_70iHqPhyK0FEzHRDyX7bSVni4lZiLWVqD2Wq8yg?key=1GDaguhNYDgftMhdZ7Y7yw)

   **Figure 2.** Heat map generated using a diverging color palette and clustered by columns


      ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdnih4QBgsiJETAw-dsZYDAuPg_4tn6plx1Gkr0gTTMxpspiDgyh0oQadlYe9T-i2nHI0AbdIj73998xV2dVbmDbb0FsOnCffkvSBLMPr8LW34eZ20C6ccGhtE6Is2fYCXtTsXhbvmdpVFcYswwhDdpD2HM?key=1GDaguhNYDgftMhdZ7Y7yw)

   **Figure 3.** Heat map generated using a diverging color palette and clustered by rows

      ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfFTnd9pUp7-6R5Dp_T6N9VmESmlBUNFRCNg_ZU2Lfxyu5LRWyBETtnn16oZGuVmHFg_p6Ch4DPxiL8VmC9JjnypY-TYsizJ9h1TAaEZ5iMNUFbNntUDYzG5ARVJGxuM4Mze3mg80bIWGx0Ew_FKRie7w?key=1GDaguhNYDgftMhdZ7Y7yw)

   **Figure 4.** Heat map generated using a diverging color palette and clustered by both columns and rows

   3.2. **Differentially expressed genes:**

   A subset of 10 overexpressed and 2 underexpressed genes were obtained.

   |          |                 |                        |                |
   | :------: | :-------------: | :--------------------: | :------------: |
   | **S.NO** |   **GENE ID**   |        **GENE**        | **EXPRESSION** |
   |     1    | ENSG00000241945 |          PWP2          | Down regulated |
   |     2    | ENSG00000279104 |         Data NA        | Down regulated |
   |     3    | ENSG00000243955 |          GSTA1         |  Up regulated  |
   |     4    | ENSG00000095917 |          TPSD1         |  Up regulated  |
   |     5    | ENSG00000231107 | Non-coding(LINC01508)  |  Up regulated  |
   |     6    | ENSG00000254092 |   Non-coding (lncRNA)  |  Up regulated  |
   |     7    | ENSG00000172236 |         TPSAB1         |  Up regulated  |
   |     8    | ENSG00000197253 |          TPSB2         |  Up regulated  |
   |     9    | ENSG00000172116 |          CD8B          |  Up regulated  |
   |    10    | ENSG00000162598 |         C1orf87        |  Up regulated  |
   |    11    | ENSG00000256193 |   Non-coding (lncRNA   |  Up regulated  |
   |    12    | ENSG00000160183 |         TMPRSS3        |  Up regulated  |

    **Table 1.** Differentially expressed genes in group 1 with respect to group 2


   3.3. **Functional enrichment analysis using ShinyGo**

   The differentially expressed genes were found to be involved in six biological pathways.

   |                                    |                                                        |                     |
   | :--------------------------------: | :----------------------------------------------------: | :-----------------: |
   |              **GENES**             |                       **PATHWAY**                      | **FOLD ENRICHMENT** |
   | TPSD1, TPSAB1,  TPSB2, and TMPRSS3 |                 Proteolysis GO:0006508                 |         5.8         |
   |                GSTA1               |  Glutathione derivative metabolic process GO:1901685   |         715         |
   |                GSTA1               | Glutathione derivative biosynthetic process GO:1901687 |         715         |
   |                 PWP2               |      Ribosomal small subunit assembly GO:0000028       |        168.2        |
   |               TMPRSS3              |       Cellular sodium ion homeostasisGO:0006883        |        124.4        |
   |                GSTA1               |       Linoleic acid metabolic process GO:0043651       |         130         |

    **Table 2.** Biological Process Pathways that were identified using ShinyGo

   ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXf3Cx4osBc7ueqwtcSHO8nC4yMGhlHdW6GDPOJ89QNNGJQJXTvi_rTp4UYCZSzx3ZRDvuusd7xEyymxv-F4Xo9yAF76NOnurbvPuwlGorOqHuWhiWygNCXY8ZxybGud98RwJr26wZIU_VlsbDJEnHnOFb4?key=1GDaguhNYDgftMhdZ7Y7yw)

    **Figure 5.** Bubble plot of function enrichment analysis

4. **DISCUSSION**

   4.1. **Importance of Colour Palette in Heatmaps**

      Choosing the right colour palette in heatmaps enhances data clarity, highlights key patterns, aids trend identification, improves accessibility for colour-blind users, and reduces cognitive load for more informative visualisations.

   4.2. **Top Three Enriched Pathways**

      4.2.1. **Glutathione derivative metabolic process**

      Pathways associated with glutathione metabolism are crucial for maintaining cellular antioxidant defense and redox balance and involve genes like the GPX family, ALDH, and SOD. These genes have been implicated in glioblastoma and other cancers (Banu et al., 2024).

      4.2.2. **Glutathione derivative biosynthetic process** 

      Focuses on the synthesis of glutathione derivatives, ensuring adequate cellular supply to counter oxidative stress and chemical challenges. (Lavoro et al., 2023) and G6PD, a protein involved in reducing oxidised glutathione, have been significantly associated with glioblastoma growth (Sun et al., 2021).

      4.3.3. **Proteolysis process**

      Regulates protein turnover, removing damaged proteins and controlling cellular functions. They include genes encoding ubiquitin proteins, caspases, and lysosomal proteins. They have associations to cancer, including glioblastoma (Auzmendi-Iriarte et al., 2022).

5. **REFERENCES** 

   Auzmendi-Iriarte, J., Otaegi-Ugartemendia, M., Carrasco-Garcia, E., Azkargorta, M., Diaz, A., Saenz-Antoñanzas, A., Andermatten, J. A., Garcia-Puga, M., Garcia, I., Elua-Pinin, A., Ruiz, I., Sampron, N., Elortza, F., Cuervo, A. M., & Matheu, A. (2022). Chaperone-Mediated Autophagy Controls Proteomic and Transcriptomic Pathways to Maintain Glioma Stem Cell Activity. _Cancer Research_, _82_(7), 1283–1297. <https://doi.org/10.1158/0008-5472.CAN-21-2161>

   Banu, M. A., Dovas, A., Argenziano, M. G., Zhao, W., Sperring, C. P., Cuervo Grajal, H., Liu, Z., Higgins, D. M., Amini, M., Pereira, B., Ye, L. F., Mahajan, A., Humala, N., Furnari, J. L., Upadhyayula, P. S., Zandkarimi, F., Nguyen, T. T., Teasley, D., Wu, P. B., … Canoll, P. (2024). A cell state-specific metabolic vulnerability to GPX4-dependent ferroptosis in glioblastoma. _The EMBO Journal_. <https://doi.org/10.1038/s44318-024-00176-4>

   Clough, E., & Barrett, T. (2016). _The Gene Expression Omnibus Database_ (pp. 93–110). <https://doi.org/10.1007/978-1-4939-3578-9_5>

   Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal patterns and correlations in multidimensional genomic data. _Bioinformatics_, _32_(18), 2847–2849. <https://doi.org/10.1093/bioinformatics/btw313>

   Lavoro, A., Falzone, L., Tomasello, B., Conti, G. N., Libra, M., & Candido, S. (2023). In silico analysis of the solute carrier (SLC) family in cancer indicates a link among DNA methylation, metabolic adaptation, drug response, and immune reactivity. _Frontiers in Pharmacology_, _14_. <https://doi.org/10.3389/fphar.2023.1191262>

   Sun, M., Sheng, H., Wu, T., Song, J., Sun, H., Wang, Y., Wang, J., Li, Z., Zhao, H., Tan, J., Li, Y., Chen, G., Huang, Q., Zhang, Y., Lan, B., Liu, S., Shan, C., & Zhang, S. (2021). PIKE-A promotes glioblastoma growth by driving PPP flux through increasing G6PD expression mediated by phosphorylation of STAT3. _Biochemical Pharmacology_, _192_, 114736. <https://doi.org/10.1016/j.bcp.2021.114736>

   Verhaak, R. G. W., Hoadley, K. A., Purdom, E., Wang, V., Qi, Y., Wilkerson, M. D., Miller, C. R., Ding, L., Golub, T., Mesirov, J. P., Alexe, G., Lawrence, M., O’Kelly, M., Tamayo, P., Weir, B. A., Gabriel, S., Winckler, W., Gupta, S., Jakkula, L., … Hayes, D. N. (2010). Integrated Genomic Analysis Identifies Clinically Relevant Subtypes of Glioblastoma Characterized by Abnormalities in PDGFRA, IDH1, EGFR, and NF1. _Cancer Cell_, _17_(1), 98–110. <https://doi.org/10.1016/j.ccr.2009.12.020> 

<!--EndFragment-->
