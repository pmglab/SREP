Readme

Here we proposed a novel approach (named SREP) to detect prognostic biomarkers.
SREP used single-cell expression clustering to purify the microenvironment of a tumor and gene pair relative expression to improve the robustness of biomarkers. 
For glioblastoma, SREP identified macrophages as the major immune cell type of the tumor and extracted differential expressed genes in macrophages. 
These genes were then paired for a subsequent prioritization for prognosis in bulk-cell gene expression data. 
In the present study, SREP identified 13 macrophage relative expression gene pairs (MREP) as prognostic biomarkers, and their accuracy and robustness were successfully validated in 3 independent multi-center datasets. 
For the first time, our study showed that usage of key genes in major immune cell types and relative gene expression substantially enhanced the accuracy and robustness of prognostic biomarkers of GBM. 
This new strategy may apply to other cancers.

We proposed a single-cell-expression-guided gene-pair prioritization workflow for detecting prognostic markers of GBM. 
It consisted of four analysis components. 

1.First, major cell types of GBM were identified using single-cell expression data of the tumor-core and tumor-peripheral tissues.
  The underlying assumption of this analysis was that the major cell types with prior knowledge-based might be more relevant to GBM and the treatment.
  The single-cells were annotated by singleR tools and revised by cellMarker database. 
  
2.Second, the maker genes list was further narrowed down by the differential expression analysis between tumor-core and tumor-peripheral tissues.
  We adopted the limma package to examine significant differentially expressed marker genes. 
  
3.Third, all retained marker genes were primarily extracted from conventional gene expression profiles and formed a gene-pair list. 
  Assume a gene pair had expression value A and B in the GBM tumor of a patient at the two genes, respectively. 
  The patient with A>B value was assigned a simplified relative expression value 1 at the gene pair; otherwise, the gene-pair value was assigned 0. 
  
4.After rigorous filtration, the UniCox was adopted to screen the gene pairs significantly associated with the prognosis of the GBM patients. 
  Finally, the significant gene pairs were further selected by LASSO-Cox (glmnet, version 4.1) and subsequently entered a multivariable Cox regression for prognostic prediction. 
  Each patient would obtain a risk score by the Cox regression, which was further used for survival prediction. 