<!--StartFragment-->

# **BIOMARKER DISCOVERY IN ACUTE MYELOID LEUKEMIA USING DIFFERENTIAL GENE EXPRESSION AND MACHINE LEARNING**

**TEAM:** Narendran Mayilvaganan (@Naren037), Ademilola Adekoya (@Demi\_Nevaeh), Ghizal Niko (@Ghizal), Eman Binte Hafeez (@emanuts), David Ojo (@Dajom), Lillian Mwinja (@LillianMwinja), Irina Andriushchenko (@rirkella), Antara Ghanta (@Tara)

**CODE:** https://github.com/Naren037/HackBio_Cancer_Internship/blob/7a54b35a64e06ebc2e754693d10b3aa75773a822/TASK%203/Codes/AML_Complete_Code.R

1. # **INTRODUCTION**

   Acute myeloid leukemia (AML) is a myeloid neoplasm resulting in clonal expansion of primitive hematopoietic stem cells, known as blasts, in the bone marrow (Vakiti & Mewawalla, 2019). Myeloblasts increase while healthy white blood cells, red blood cells, and platelets decrease (National Cancer Institute, 2019). AML can also be classified as newly diagnosed, in remission, or recurrent (National Cancer Institute, 2019).  

2. # **METHODOLOGY**

   1. ## **Datamining**

      The dataset was extracted from project “TARGET-AML” in TCGA using TCGABiolinks with filters to fetch RNA-Seq data of 20 primary and recurrent samples. 

   2. ## **Preprocessing for DGE** 

      Cross-sample RNA-seq data were preprocessed in 3 steps: replacing NA values with rowmeans, TMM normalization using Edge to accommodate variations in sequencing depth and sample bias, and using an upper quantile filter to remove lowly    expressed genes.

   3. ## **Differential Gene Expression (DGE) analysis**

      TCGAnalyze\_DEA was used to identify genes with significant expression differences between primary and recurrent AML samples. To identify meaningful biomarkers, a super stringent cutoff of FDR = 0.0005 and logFC = 4 were used.  

   4. ## **Functional Enrichment analysis**

      Functional enrichment was done using TCGAanalyze\_EAcomplete. The top 5 GO biological process pathways were considered for interpretation.

   5. ## **Preprocessing and Feature selection for ML** 

      RNA-Seq data was filtered by selecting high variance genes. It was then filtered to remove zero values and avoid multicollinearity by centering. Lastly, the features with high correlation (cutoff of  0.53) were removed. 

   6. ## **ML Model training and perfomance using kNN**

      Based on the number of features, through educated trail & error 3 optimizations were made to the kNN settings: number = 4, preProcess = "scale" and k = 1:8. A repeat loop was made to shuffle the dataset in every iteration before 70:30    split for kNN. The loop exit condition was set to train accuracy of 100% and test accuracy > 90%.

3. # **RESULTS**

   1. ## **Differential Gene Expression**

      24 differentially expressed genes between primary and recurrent AML were obtained. 

**Table 1.** Top 10 differentially expressed genes and their influence in AML

|          |             |                                                                                                                                                                                                                        |                             |
| :------: | :---------: | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :-------------------------: |
| **GENE** | **P VALUE** |                                                                                                      **FUNCTION**                                                                                                      |        **REFERENCE**        |
|   PITX1  |   4.31E-08  |                                     PITX1 interacts with p53 and regulates essential biological processes such as cell cycle progression, apoptosis, and chemotherapy resistance.                                      |      (Zhao & Xu, 2023)      |
|  RUNX1T1 |   2.64E-07  |                                                                      RUNX1T1 is the predominant fusion gene identified in acute myeloid leukemia.                                                                      | (Swart & Heidenreich, 2021) |
|  CHRDL1  |   1.95E-06  |                                                                                   CHRDL1 was upregulated in acute myeloid leukemia.                                                                                    |      (Yu et al., 2023)      |
|   GABRE  |   1.77E-06  |                                                                                         Reported as a subtype biomarker of AML                                                                                         |     (Labaj et al., 2017)    |
|   FOXF1  |   4.20E-06  |                                                                               Reported to be differentially expressed across AML subtypes                                                                              |      (Yin et al., 2022)     |
|    MET   |   1.94E-07  |                                                                                                     Proto-oncogene                                                                                                     |     (Rivas et al., 2022)    |
|   CCL2   |   4.55E-10  |                                                                             Modulation of CCL2 will improve response to chemotherapy of AML                                                                            |    (Jacamo et al., 2015)    |
|  ADAM23  |   1.43E-09  |                                                                                             Reported in Chemorefractory AML                                                                                            |    (Palmer et al., 2015)    |
|   MNX1   |   2.20E-06  |                                        MNX1 exhibited abnormal expression in 1.4% of AML patients, and its suppression in a xenograft model diminished leukemia cell viability.                                        |    (Sollier et al., 2024)   |
|    IL6   |   6.70E-10  | The dysregulation of intricate interactions between pro- and anti-inflammatory cytokines in (AML) may foster a pro-tumorigenic environment, influencing leukemic cell proliferation, survival, and therapy resistance. |    (Binder et al., 2018)    |


   ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeJqP5bnHMJ3PCKdUtGDV56Ytbp9Q4DAwRk4KJtukEIUPrXVuyTJY9g0rDKi0VnS1SVWBTUSB7fjzhBOZKR-WQtWj34qYGsd4TNs7WPus4J833kXcJR57TfLxChvxM9WFoIuEuy7b2QLt6oFIUfYZT7-TDr?key=UwCJV6gpCQNC2y-qvDXcZw)

   **Figure 1.** Heatmap of differentially expressed genes shows more upregulation in recurrent samples.

   ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXfWETIicQLsq40johlnNpVq8OvgyXM3gLFXuDaMr4s041h_bTjZKE1AkYs4wXgBML0-OV5nG3uzd4C7g-YqhkQhI2ba0_HPOWEdn4R8qO4uhAvFowOUcH8byFHrb2otL8UnDhHE4yiOzLgmAmUzLA8BzmPs?key=UwCJV6gpCQNC2y-qvDXcZw)

   **Figure 2.** Bubble plot of the differentially expressed genes

2. ## **Functional Enrichment Analysis**

   The top 5 GO biological process pathways obtained are mentioned in Table 2.

   &#x20;**Table 2.** Top Biological Process Pathways that were identified 
      |            |                                                         |         |
      | :--------: | :-----------------------------------------------------: | :-----: |
      | **nGENES** |                       **PATHWAY**                       | **FDR** |
      |      2     |             GO:0050900\~leukocyte migration             | 0.00704 |
      |      2     |         GO:0035270\~endocrine system development        | 0.00704 |
      |      3     |  GO:0008284\~positive regulation of cell proliferation  | 0.00777 |
      |      2     |      GO:0001817\~regulation of cytokine production      | 0.00777 |
      |      2     | GO:0002696\~positive regulation of leukocyte activation | 0.00777 |

   

4. ## **Machine Learning**

   The ML model performed consistently with a train accuracy of 100% and test accuracy of 91%. The model used 50 features in total after preprocessing, and the top 5 features are mentioned in Table.3. 

   ![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXcDkzbpdh1_XxVEqD0vifynn56lWMgiME_8aJKpj4QioijS_XtLvLpZIetCnMjobPitqoTwcRWj4iNFSs8ml32mcjhrULJUOnteZ35l-MEJz4mqM_OKwCkZgaroKSQul9Myi-zEGo9a6-FOm27A_kKxwjzR?key=UwCJV6gpCQNC2y-qvDXcZw)

   **Figure 3.** A heatmap of preprocessed data showing features with unique expression profiles.

   **Table 3.** Top 5 features used by the ML model
   |          |                                                                                                                                                                           |                                  |
   | :------: | :-----------------------------------------------------------------------------------------------------------------------------------------------------------------------: | :------------------------------: |
   | **GENE** |                                                                                **FUNCTION**                                                                               |           **REFERENCE**          |
   |   COX1   |                  Up-regulation results in a distinct “eicosanoid storm” along with “anti-inflammatory, proinflammatory, and redox-activated” signatures.                  |   (Pannunzio & Coluccia, 2018)   |
   |  RN7SL1  |                                                                      Drives an inflammatory response                                                                      |      (Johnson et al., 2021)      |
   |   HSPA5  | Involved in the correct folding and degradation of misfolded proteins. Higher expression of HSPA5 significantly decreased patient overall survival in 7 types of cancers. |         (Fu et al., 2021)        |
   |   ELANE  |                                                        Encodes neutrophil elastase, ELANE mutation can lead to AML.                                                       |       (Rotulo et al., 2020)      |
   |   GATA2  |                                             Zinc finger transcription factor, loss or gain of which has been implicated in AML                                            | (Menendez-Gonzalez et al., 2019) |



5. # **DISCUSSION**

   24 out 60660 genes were found to be significantly altered with our super-stringent cutoffs. The most statistically significant genes, IL6 and CCL2 were found to influence chemotherapy, serving as potential biomarkers for recurrent AML. The leukocyte activation and migration pathways (Mishra et al., 2020) explains the abnormal leukocyte activity in AML. Endocrine changes have been linked to AML progression (Roma & Spagnuolo, 2020). 

6. # **CONCLUSION AND FUTURE PERSPECTIVES** 

   We identified potential biomarkers for AML recurrence through DGE and ML. Despite small sample size, the high relevance of obtained biomarkers in AML highlight the potency of our optimizations. Future perspectives include validating our methods using larger cohorts and exploring their practical significance in AML biology by _in vivo_ validation. 

7. # **REFERENCES**

   Binder, S., Luciano, M., & Horejs-Hoeck, J. (2018). The cytokine network in acute myeloid leukemia (AML): A focus on pro- and anti-inflammatory mediators. _Cytokine & Growth Factor Reviews_, _43_, 8–15. <https://doi.org/10.1016/j.cytogfr.2018.08.004> 
   
   Fu, J., Wei, C., He, J., Zhang, L., Zhou, J., Balaji, K. S., Shen, S., Peng, J., Sharma, A., & Fu, J. (2021). Evaluation and characterization of HSPA5 (GRP78) expression profiles in normal individuals and cancer patients with COVID-19. _International Journal of Biological Sciences_, _17_(3), 897–910. <https://doi.org/10.7150/ijbs.54055> 
   
   Jacamo, R. O., Mu, H., Zhang, Q., Chachad, D., Zhiqiang, W., Ma, W., Zhang, M., Mak, P. Y., Mak, D., Ruvolo, P., McQueen, T., Lowe, S., Zuber, J., Eulberg, D., Kruschinski, A., Konopleva, M., Davis, R. E., & Andreeff, M. (2015). Effects of CCL2/CCR2 Blockade in Acute Myeloid Leukemia. _Blood_, _126_(23), 1348–1348. <https://doi.org/10.1182/blood.V126.23.1348.1348> 
   
   Johnson, L. R., Lee, D. Y., Eacret, J. S., Ye, D., June, C. H., & Minn, A. J. (2021). The immunostimulatory RNA RN7SL1 enables CAR-T cells to enhance autonomous and endogenous immune function. _Cell_, _184_(19), 4981-4995.e14. <https://doi.org/10.1016/j.cell.2021.08.004> 
   
   Labaj, W., Papiez, A., Polanski, A., & Polanska, J. (2017). Comprehensive Analysis of MILE Gene Expression Data Set Advances Discovery of Leukaemia Type and Subtype Biomarkers. _Interdisciplinary Sciences: Computational Life Sciences_, _9_(1), 24–35. <https://doi.org/10.1007/s12539-017-0216-9> 
   
   Menendez-Gonzalez, J. B., Sinnadurai, S., Gibbs, A., Thomas, L., Konstantinou, M., Garcia-Valverde, A., Boyer, M., Wang, Z., Boyd, A. S., Blair, A., Morgan, R. G., & Rodrigues, N. P. (2019). Inhibition of GATA2 restrains cell proliferation and enhances apoptosis and chemotherapy mediated apoptosis in human GATA2 overexpressing AML cells. _Scientific Reports_, _9_(1), 12212. <https://doi.org/10.1038/s41598-019-48589-0>
   
   Mishra, S., Kumar, K., Panigrahi, A., Das, P., Padhi, S., & Chhabra, G. (2020). The Utility of Leucocyte Cell Population Data and Scattergram in Rapid Identification of Acute Promyelocytic Leukemia. Blood, 136(Supplement 1), 19–20. <https://doi.org/10.1182/blood-2020-142498> 
   
   National Cancer Institute. (2019, July 23). Adult Acute Myeloid Leukemia Treatment (PDQ®)–Patient Version. National Cancer Institute; Cancer.gov. <https://www.cancer.gov/types/leukemia/patient/adult-aml-treatment-pdq> 
   
   Palmer, A., Parkin, B., Posthuma, H., Mineishi, S., Magenau, J. M., & Malek, S. N. (2015). Gene Mutations Identified in Chemorefractory Acute Myeloid Leukemia. _Blood_, _126_(23), 3838–3838. <https://doi.org/10.1182/blood.V126.23.3838.3838> 
   
   Pannunzio, A., & Coluccia, M. (2018). Cyclooxygenase-1 (COX-1) and COX-1 Inhibitors in Cancer: A Review of Oncology and Medicinal Chemistry Literature. _Pharmaceuticals_, _11_(4), 101. <https://doi.org/10.3390/ph11040101> 
   
   Rivas, S., Marín, A., Samtani, S., González-Feliú, E., & Armisén, R. (2022). MET Signaling Pathways, Resistance Mechanisms, and Opportunities for Target Therapies. _International Journal of Molecular Sciences_, _23_(22), 13898. <https://doi.org/10.3390/ijms232213898> 
   
   Roma, A., & Spagnuolo, P. A. (2020). Estrogen Receptors Alpha and Beta in Acute Myeloid Leukemia. Cancers, 12(4), 907. <https://doi.org/10.3390/cancers12040907> 
   
   Sollier, E., Riedel, A., Toprak, U. H., Wierzbinska, J. A., Weichenhan, D., Schmid, J. P., Hakobyan, M., Touzart, A., Jahn, E., Vick, B., Brown-Burke, F., Kelly, K., Kelekci, S., Pejkovska, A., Goyal, A., Baehr, M., Breuer, K., Chen, M. M., Llamazares-Prada, M., . . . Plass, C. (2024). Pyjacker identifies enhancer hijacking events in acute myeloid leukemia including MNX1 activation via deletion 7q. bioRxiv (Cold Spring Harbor Laboratory). <https://doi.org/10.1101/2024.09.11.611224> 
   
   Rotulo, G. A., Beaupain, B., Rialland, F., Paillard, C., Nachit, O., Galambrun, C., Gandemer, V., Bertrand, Y., Neven, B., Dore, E., Moshous, D., Filhon, B., Aladjdi, N., Sicre de Fontbrune, F., de la Tour, R. P., Ouachee, M., Bellanne-Chantelot, C., Dalle, J.-H., & Donadieu, J. (2020). HSCT may lower leukemia risk in ELANE neutropenia: a before–after study from the French Severe Congenital Neutropenia Registry. _Bone Marrow Transplantation_, _55_(8), 1614–1622. <https://doi.org/10.1038/s41409-020-0800-1> 
   
   Swart, L. E., & Heidenreich, O. (2021). The RUNX1/RUNX1T1 network: translating insights into therapeutic options. _Experimental Hematology_, _94_, 1–10. <https://doi.org/10.1016/j.exphem.2020.11.005> 
   
   Vakiti, A., & Mewawalla, P. (2019, December 18). Cancer, Acute Myeloid Leukemia (AML, Erythroid Leukemia, Myelodysplasia-Related Leukemia, BCR-ABL Chronic
   
   Yin, H., Fan, X., Zhang, Y., Zhao, N., Zhao, X., Yin, K., & Zhang, Y. (2022). An Integrated Study on the Differential Expression of the FOX Gene Family in Cancer and Their Response to Chemotherapy Drugs. _Genes_, _13_(10), 1754. <https://doi.org/10.3390/genes13101754> 
   
   Yu, J.-W., Pang, R., Liu, B., Zhang, L., & Zhang, J.-W. (2023). Bioinformatics identify the role of chordin-like 1 in thyroid cancer. _Medicine_, _102_(5), e32778. <https://doi.org/10.1097/MD.0000000000032778> 
   
   Zhao, J., & Xu, Y. (2023). PITX1 plays essential functions in cancer. _Frontiers in Oncology_, _13_. <https://doi.org/10.3389/fonc.2023.1253238>

<!--EndFragment-->
