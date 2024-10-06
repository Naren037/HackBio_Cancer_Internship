<!--StartFragment-->

# **“BIOMARKER DISCOVERY IN LOW GRADE GLIOBLASTOMA USING DIFFERENTIAL GENE EXPRESSION AND UNSUPERVIZED MACHINE LEARNING”**

**TEAM:** Narendran Mayilvaganan (@Naren037), Ademilola Adekoya (@Demi\_Nevaeh), Ghizal Niko (@Ghizal), Eman Binte Hafeez (@emanuts), David Ojo (@Dajom), Lillian Mwinja (@LillianMwinja), Irina Andriushchenko (@rirkella), Antara Ghanta (@Tara)

**CODE:** <https://github.com/Naren037/HackBio_Cancer_Internship/blob/740966569cf42bdd40078827a987517144062de8/TASK%203/Codes/AML_Complete_Code.R> 

1. # **INTRODUCTION**<a id="h.yv49cqze5iy3"></a>

   A glioma is a tumor of the central nervous system that originates from glial cells and is the most common brain tumor in adults (Wesseling & Capper, 2018). Until recently, the only factor       used to classify gliomas was their morphological characteristics; however, IDH mutations are now considered to be crucial in this process. This classification has a prognostic value as           patients with IDH-mut tumors had higher survival rates than those with IDH-wt tumors (Hua et al., 2024).

2. # **METHODOLOGY**<a id="h.9z72xc2z3mv2"></a>

      1. ## **Datamining**<a id="h.4n3jiudp7xpj"></a>

         The dataset was extracted from project “TARGET-LGG” in TCGA using TCGABiolinks with filters to fetch RNA-Seq data of n=534 LGG samples.

      2. ## **Preprocessing for DGE**<a id="h.f9lrawf6ktzt"></a>

         Cross-sample RNA-seq data were preprocessed in 3 steps: replacing NA values with rowmeans, TMM normalization using EdgeR to accommodate variations in sequencing depth and sample bias,            and using an upper quantile filter to remove lowly expressed genes.

      3. ## **Differential Gene Expression (DGE) analysis**<a id="h.m74rnaszgj4"></a>

         TCGAnalyze\_DEA was used to identify genes with significant expression differences between WT and Mutant-IDH samples. To identify meaningful biomarkers, a super stringent cutoff of FDR           = 0.0005 and logFC = 4 were used.  

      4. ## **Functional Enrichment analysis**<a id="h.iseit6flhn21"></a>

         Functional enrichment was done using TCGAanalyze\_EAcomplete. The top 5 GO biological process pathways were considered for interpretation.

      5. ## **Preprocessing and Feature selection for ML** <a id="h.gyjoigaekn3t"></a>

         RNA-Seq data was normalized using DESeq2 and filtered by selecting high variance genes. It was then filtered to remove zero values and avoid multicollinearity by centering. Lastly, the           features with high correlation (cutoff of  0.5) were removed. 

      6. ## **Unsupervised Machine Learning**<a id="h.p2r45tt8njbp"></a>

         Dimensionality reduction was done using PCA with filtered features. The PCA plot was then combined with the clustering data obtained from unsupervised ML using k-Means. The number of             clusters were double-checked using Hierarchical clustering. Finally the IDH subtypes of each sample in the clusters were retrieved and studied. 

3. # **RESULTS**<a id="h.1bwiz9ej2jbz"></a>

   1. ## **DGE**<a id="h.246f2h4ntzsv"></a>

      34 differentially expressed genes between WT and recurrent IDH-Mutant were obtained.

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcsgF38gjWGpLscBM1Fy0yfUYAlZjhkDDZFwqNdq5zDDJ-Mta_5532HUk0O7or6K0H2J5Y3Ugwqhul1wcUACObOpL_ccPcE3mo6s3I-CVygrjnoDmXTTh5J5EaWB5VGjfmFu12P0hjjmNiol5oRjhYKofQ?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 1**. Heatmap of DEGs show similar cluster patterns for IDH Mutant and WT samples respectively

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcuz0IL0PFq6Lmp5kfBD5EVHpYu2CP455Uin2WooqnA-beGg9hvB5t6iT5P9SyrXAOXZkl-uXcMfr_QEhdfhOb585ikM78dqMWi5T-5N12-gq7oUTY5B5mVzaRQA1rZMkdZABFzdh7s_aISCTTJErMwUQI?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 2.** Bubble plot of the differentially expressed genes

   2. ## **Functional enrichment**<a id="h.qhy625fyw8fy"></a>

|            |                                                   |         |
| :--------: | :-----------------------------------------------: | :-----: |
| **nGENES** |                    **PATHWAY**                    | **FDR** |
|      2     |     GO:0048562\~embryonic organ morphogenesis     | 0.00800 |
|      2     |      GO:0048568\~embryonic organ developmentt     | 0.00800 |
|      1     |    GO:0048617\~embryonic foregut morphogenesis    | 0.00800 |
|      2     |     GO:0007389\~pattern specification process     | 0.00800 |
|      2     | GO:0009792\~embryonic development ending in birth | 0.00800 |

&#x20;**Table 1.** Top 5 GO Biological Processes that were identified 

   3. ## **Unsupervized ML - PCA k means and clusterring**<a id="h.t45me9ebfogm"></a>

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeiq5YSfZUqTsPSOiNk89VkiwtbuLO5xmgAXdnLX8q3HShVW8yE1x7XjD6iUHhM_kbFbPiXzgudANlENOhC-Se3qAR10fJ39P3K8L66-yW-rwQNCCeGSmULBgsZSccMCiLqBXXcKxR83iCf9C1MC9boqf_o?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 3.** PCA visualization of the preprocessed dataset 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeJIYixWNKiSnkkWHELirPDujZBgEW-AARQgG9177gHNx_o8mwLr6rq4zZIBSZHV0qH2JH64KKu0VDKD0kUoa-Q94uiqz2JGER5D7yT7aPX2dUw43Y-GLrWFINQ9wiRXvLsex9OuvlY7VTGhgbU4CYFI7Y?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 4.** PCA visualization after unsupervised ML using k-Means methods show 3 distinct cluster groups

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdEv29ItM15WctDvVUNONZVNOA-1uUQEVnlRXXogl9EXu6Qmpc88XrXsfjkOTSS6HuN17XcnU7bD7MjcJMJojNCQQxVbLdV7QYOYSbi-4dSfX_Duflv3xg3wtNHBkEu4j2Uapw0UUFlp-nHxL0LTcrJwu9X?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 5.** Unsupervised ML using Hierarchical clustering was done to further verify the number of major cluster groups

**Table 2.** No.of samples predicted in each cluster for IDH mutations  

|                      |               |               |               |
| -------------------- | ------------- | ------------- | ------------- |
| **Mutation Subtype** | **Cluster 1** | **Cluster 2** | **Cluster 3** |
| IDHmut-codel         | 0             | 160           | 9             |
| IDHmut-non-codel     | 234           | 2             | 14            |
| IDHwt                | 0             | 0             | 94            |

**Table 3.** No.of samples predicted in each cluster for ATRX mutations 

|                      |               |               |               |
| -------------------- | ------------- | ------------- | ------------- |
| **Mutation Subtype** | **Cluster 1** | **Cluster 2** | **Cluster 3** |
| Mutant               | 173           | 2             | 13            |
| WT                   | 61            | 160           | 104           |

**Table 4.** No.of samples predicted in each cluster for TERT mutations

|                      |               |               |               |
| -------------------- | ------------- | ------------- | ------------- |
| **Mutation Subtype** | **Cluster 1** | **Cluster 2** | **Cluster 3** |
| Expressed            | 10            | 148           | 66            |
| Not expressed        | 226           | 14            | 49            |

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc3N64LlxX0P35k7UAvdv6ZLoUPTujudMyhVERzZ13xzl1_cCUlIwbdXWgT6dz7WUBUV5ZoHGfNNyybDeEui-4WlK6mSoUHSSIR23Ux_WFc1QidHwsk2V4E6OAniYD14zgPm5-zZ4fC8DDLJetNHXMqrwq7?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 6.** Barplot of IDH Mutations in each cluster

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc5YtT1jz9o2MRofV-19JUVfUbUpSDGSF020RNODulp1lEIdyDWIwaLZK5xUkTUSyB-KCg459Zc14h7NWXtY2GesXPjSkNCSQ7FFzzkZAXxD2kZcGsy4mYHG2fooRAmXHODrZqrW3HilDauWph72aXYMJs?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 7.** Barplot of ATRX Mutations in each cluster

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXdcd9mKvQqrGQwyrWa0a8suArt8u3nbY2tb8WWlOCbD2zkPQxJ5jSLPn-3hCOyNJSMLFgkf4CdrpF9pFdqQB9OnteV6cqBm5UwmJkpa5b5NuSHsfSaJzmQ7uoO_gyBmwthzEqc5vQqsvuX8GrolMROeCloT?key=83q9Nzh7nBQo6FjuiFK5mg)

**Figure 8.** Barplot of TERT Mutations in each cluster

4. # **DISCUSSION**<a id="h.4ceafm4pi6lh"></a>

   The unsupervised ML clustered the gene expression data into 3 groups with high accuracy. The clusters were found to be majorly distributed as, 

   Cluster 1- IDHwt & ATRXwt

   Cluster 2 - IDHmut non-codel & TERT Not expressed 

   Cluster 3 - IDHmut codel & ATRXwt & TERT Expressed

   In comparison, the Ceccarelli’s study (Ceccarelli  et al., 2016) identified six methylation groups and four RNA expression groups associated with IDH status. Our ML findings are similar to       the results obtained in Ceccarelli’s study : TERT is mostly not expressed when there is also a IDH mut non-codel and expressed in case of IDH mut codel. 34 out 60660 genes were found to be       significantly altered with our super-stringent cutoffs. The most statistically significant genes, MEOX2, FMOD, LTF, RARRES2, DMRTA2 and MAP1LC3C (Guo et al., 2021; Mondal et al., 2017;           Schönrock et al., 2022) were found to be associated with LGG. Interestingly, few non-coding genes were reported with with >10 Foldchange. The Functional Enrichment majorly revealed               embryogenesis and organ development GO-processes. We conclude that more research needs to be done to unravel the role of non-coding genes in LGG and that IDH mutations may affect organ           development, increaing susceptibility to LGGs.

5. # **REFERENCES**<a id="h.breuynwez23"></a>

   Ceccarelli, M., Barthel, F. P., Malta, T. M., Sabedot, T. S., Salama, S. R., Murray, B. A., Morozova, O., Newton, Y., Radenbaugh, A., Pagnotta, S. M., Anjum, S., Wang, J., Manyam, G.,            Zoppoli, P., Ling, S., Rao, A. A., Grifford, M., Cherniack, A. D., Zhang, H., … Zmuda, E. (2016). Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in         Diffuse Glioma. _Cell_, _164_(3), 550–563. <https://doi.org/10.1016/j.cell.2015.12.028>

   Guo, J.-C., Wei, Q.-S., Dong, L., Fang, S.-S., Li, F., & Zhao, Y. (2021). Prognostic Value of an Autophagy-Related Five-Gene Signature for Lower-Grade Glioma Patients. _Frontiers in              Oncology_, _11_. <https://doi.org/10.3389/fonc.2021.644443>

   Hua, W., Zhang, W., Brown, H., Wu, J., Fang, X., Shahi, M., Chen, R., Zhang, H., Jiao, B., Wang, N., Xu, H., Fu, M., Wang, X., Zhang, J., Zhang, X., Wang, Q., Zhu, W., Ye, D., Garcia, D. M.,    … Quinones-Hinojosa, A. (2024). Rapid detection of IDH mutations in gliomas by intraoperative mass spectrometry. _Proceedings of the National Academy of Sciences_, _121_(23).                   
   <https://doi.org/10.1073/pnas.2318843121>

   Mondal, B., Patil, V., Shwetha, S. D., Sravani, K., Hegde, A. S., Arivazhagan, A., Santosh, V., Kanduri, M., & Somasundaram, K. (2017). Integrative functional genomic analysis identifies         epigenetically regulated fibromodulin as an essential gene for glioma cell migration. _Oncogene_, _36_(1), 71–83. <https://doi.org/10.1038/onc.2016.176>

   Schönrock, A., Heinzelmann, E., Steffl, B., Demirdizen, E., Narayanan, A., Krunic, D., Bähr, M., Park, J.-W., Schmidt, C., Özduman, K., Pamir, M. N., Wick, W., Bestvater, F., Weichenhan, D.,     Plass, C., Taranda, J., Mall, M., & Turcan, Ş. (2022). _MEOX2_ homeobox gene promotes growth of malignant gliomas. _Neuro-Oncology_, _24_(11), 1911–1924.                      
   <https://doi.org/10.1093/neuonc/noac110>

   Wesseling, P., & Capper, D. (2018). \<scp>WHO\</scp> 2016 Classification of gliomas. _Neuropathology and Applied Neurobiology_, _44_(2), 139–150. <https://doi.org/10.1111/nan.12432> 

<!--EndFragment-->
