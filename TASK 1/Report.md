<!--StartFragment-->


# **MULTI-OMIC MACHINE LEARNING PREDICTOR OF BREAST CANCER THERAPY RESPONSE**

**Stephen-John Sammut**, Mireia Crispin-Ortuzar, Suet-Feung Chin, Elena Provenzano, Helen A. Bardwell, Wenxin Ma, Wei Cope, Ali Dariush, Sarah-Jane Dawson, Jean E. Abraham, Janet Dunn, Louise Hiller, Jeremy Thomas, David A. Cameron, John M. S. Bartlett, Larry Hayward, Paul D. Pharoah, Florian Markowetz, Oscar M. Rueda, Helena M. Earl & Carlos Caldas.

**DOI**: <https://doi.org/10.1038/s41586-021-04278-5> 

**LITERATURE REVIEW**

**Authors:** Eman Binte Hafeez (@emanuts**)**, Fatma Fathi (@FatmaFathi), Ghizal Niko (@Ghizal) and Narendran Mayilvaganan (@Naren037).


# **1. INTRODUCTION**

Neoadjuvant treatment, involving chemotherapy with or without targeted therapy before surgery, is used to        improve breast cancer outcomes but often has limited effectiveness (Symmans et al., 2007, 2017). Clinical features are not conclusive indicators, and other factors do influence treatment response (Marusyk et al., 2020). The study used multi-omics data from pre-therapy biopsies to build a machine-learning prediction model integrating clinical, molecular, and pathology data.


# **2. MULTI-OMIC FEATURES AND MACHINE LEARNING PREDICTION**

## **2. 1 Multi-platform profiling of tumour biopsies**

The clinical data was obtained from 180 women with early and locally advanced breast cancer undergoing neoadjuvant treatment (refer Figure 1). Pre-treatment tumour biopsies were analysed using whole-genome, whole-exome, and RNA sequencing. Chemotherapy was administered to all patients, whereas 65 patients with HER2+ tumours received further anti-HER2 gene therapy. 

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXc8oCBnpOS5kQng1kd_EGEGzO92GtbJF6N5MQXrbvN_Is_DKpT3-J_K_Z3SckKl-Kbkx9uKEAHC5_FvNitaPgrTWRiE4gLC2dbjz6bmM_yLsx8LpkZnZRHKOCXpIyrF5eOGIQjQfkw8k-WbC879f13s0O6y?key=WjCHy4SolYKCRAPkqrfFyw)

**Figure 1**. Clinical workflow and data acquisition (Sammut et al., 2022)

**Table 1**.Response assessment at surgery using the residual cancer burden (RCB) classification.

| Response (RCB)                       | Assessment  |
| ------------------------------------ | :---------: |
| Pathological Complete Response (pCR) |     26%     |
| Good response (RCB-I)                |     16%     |
| Moderate response (RCB-II)           |     40%     |
| Extensive residual disease (RCB-III) |     18%     |


## **2. 3 Genomic landscapes and other factors associated with response**

Whole-exome sequencing identified 16,134 somatic mutations, with certain major mutations that are mentioned in Table 2. 

**Table 2**. Major mutations identified

|                  |      |        |       |        |
| :--------------: | :--: | :----: | :---: | :----: |
| **Mutated gene** | TP53 | PIK3CA | GATA3 | MAP3K1 |
|  **Frequency**   |  57% |   26%  |  10%  |   8%   |

pCR tumours:

- Greater chromosomal instability, more copy number alterations, higher mutation and neoantigen burdens, and simpler clonal structures. 

- P53 mutations and specific genomic markers APOBEC and HRD.

Resistant tumours:

- Tumours lacking the HLA-I gene.

- PIK3CA mutations.


## **2.4 Tumour proliferation, TiME, and immune system attenuation**

Differential RNA expression study revealed 2,071 genes underexpressed and 2,439 overexpressed genes in pCR tumours. Gene-set enrichment analysis (Liberzon et al., 2015), showed a strong correlation between proliferation and immune activation with pCR. GGI analysis (Sotiriou et al., 2006) showed tumours that attained pCR had higher proliferation. Tumour immunological Micro Environment (TiME) was analysed using RNA expression deconvolution, showing enriched innate and adaptive immunity in pCR. Whereas, resistant tumours have increased mast cells. Active TiME didn’t guarantee pCR due to an increase in inhibitory natural killer CD56 cells and regulatory T-cells, and increased exclusion in resistant tumours.


## **2.5 Integration of multi-omics features with Machine Learning** 

The ML models were developed by using only clinical pathology data and subsequently expanding the number of features. The pipeline was created using three parallel algorithms to analyse the parameters in each model and was tested with a fivefold cross-validation scheme.

****![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXehbEadHaE3bJKh7Is--qn4-ma3ZUzJHTqDLISAp_7NxfhPQ8bIlEde4o8IrFQ9Tlb_pcrtJ3Oi4g-swNXM8hhxBDT5nHbwS2iRavBiImpL-7kvVzZK80rLr1h560caOeSWKy6m8-T967LV83qH3x8T62I?key=WjCHy4SolYKCRAPkqrfFyw)**                                      Figure 2**. Pipeline for prediction using ML.********


# **3. METHODOLOGY**

Blood samples from patients with primary invasive cancer were collected before starting therapy, and RCB, ER, and HER2 were assessed. DNA quantification and sequencing, and RNA extraction and transcriptomics were done (Dobin et al., 2013). A bioinformatics pipeline was made for preprocessing, clonal reconstruction, immune-related analysis, tumour classification, RNA-seq alignment, and gene counts obtained using HTseq.


# **4. LIMITATIONS**

Ambiguity due to heterogeneity due to the small sample size of 168 tumours. Whole-exome sequencing overlooks epigenetic alterations. Findings limited only to certain tumour subtypes. Mutations and treatment responses are correlational, not causal. There are immune escape mechanisms other than HLA-I. Limited variability in treatment and lack of long-term follow-up data.


# **5. DISCUSSION AND CONCLUSION**

The study demonstrates the use of ML to predict breast cancer therapy responses by combining clinical, molecular, and digital pathology to accurately predict pathological responses. Therefore, the use of data analysis and ML for personalised medicine can circumvent the inability of clinical features to be used as prediction tools in cancer.


# **6. REFERENCES**

Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: ultrafast universal RNA-seq aligner. _Bioinformatics_, _29_(1), 15–21. https\://doi.org/10.1093/bioinformatics/bts635

Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M., Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database Hallmark Gene Set Collection. _Cell Systems_, _1_(6), 417–425. https\://doi.org/10.1016/j.cels.2015.12.004

Marusyk, A., Janiszewska, M., & Polyak, K. (2020). Intratumor Heterogeneity: The Rosetta Stone of Therapy Resistance. _Cancer Cell_, _37_(4), 471–484. https\://doi.org/10.1016/j.ccell.2020.03.007

Sammut, S.-J., Crispin-Ortuzar, M., Chin, S.-F., Provenzano, E., Bardwell, H. A., Ma, W., Cope, W., Dariush, A., Dawson, S.-J., Abraham, J. E., Dunn, J., Hiller, L., Thomas, J., Cameron, D. A., Bartlett, J. M. S., Hayward, L., Pharoah, P. D., Markowetz, F., Rueda, O. M., … Caldas, C. (2022). Multi-omic machine learning predictor of breast cancer therapy response. _Nature_, _601_(7894), 623–629. https\://doi.org/10.1038/s41586-021-04278-5

Sotiriou, C., Wirapati, P., Loi, S., Harris, A., Fox, S., Smeds, J., Nordgren, H., Farmer, P., Praz, V., Haibe-Kains, B., Desmedt, C., Larsimont, D., Cardoso, F., Peterse, H., Nuyten, D., Buyse, M., Van de Vijver, M. J., Bergh, J., Piccart, M., & Delorenzi, M. (2006). Gene Expression Profiling in Breast Cancer: Understanding the Molecular Basis of Histologic Grade To Improve Prognosis. _JNCI: Journal of the National Cancer Institute_, _98_(4), 262–272. https\://doi.org/10.1093/jnci/djj052

Symmans, W. F., Peintinger, F., Hatzis, C., Rajan, R., Kuerer, H., Valero, V., Assad, L., Poniecka, A., Hennessy, B., Green, M., Buzdar, A. U., Singletary, S. E., Hortobagyi, G. N., & Pusztai, L. (2007). Measurement of Residual Breast Cancer Burden to Predict Survival After Neoadjuvant Chemotherapy. _Journal of Clinical Oncology_, _25_(28), 4414–4422. https\://doi.org/10.1200/JCO.2007.10.6823

Symmans, W. F., Wei, C., Gould, R., Yu, X., Zhang, Y., Liu, M., Walls, A., Bousamra, A., Ramineni, M., Sinn, B., Hunt, K., Buchholz, T. A., Valero, V., Buzdar, A. U., Yang, W., Brewster, A. M., Moulder, S., Pusztai, L., Hatzis, C., & Hortobagyi, G. N. (2017). Long-Term Prognostic Risk After Neoadjuvant Chemotherapy Associated With Residual Cancer Burden and Breast Cancer Subtype. _Journal of Clinical Oncology_, _35_(10), 1049–1060. <https://doi.org/10.1200/JCO.2015.63.1010>

Videolink: <https://www.linkedin.com/posts/eman-binte-hafeez-00a455283_hackbio-breastcancerresearch-multiomics-activity-7238214931987238912-74WD?utm_source=share&utm_medium=member_desktop>

\
\


<!--EndFragment-->
