# Semi-supervised-learning-for-cell-type-identification-MSc-thesis-
This repository contains the code I have written for my MSc thesis on semi-supervised learning for cell type identification using scRNA-seq data.

## Abstract
Recent research in bioinformatics has presented multiple cell type identification meth- dologies using single cell RNA sequence data (scRNA-seq). However, a consensus on which cell typing methodology consistently demonstrates superior performance remains absent. Additionally, very few studies approach cell type identification through a semi- supervised learning study, whereby the information in unlabeled data is leveraged to train an enhanced classifier. This paper presents cell annotation methodologies through self- learning and graph-based semi-supervised learning, in both raw count scRNA-seq data as well as in a latent embedding. I find that a self-learning framework enhances perfor- mance compared to a solely supervised learning classifier. Additionally, modelling on the latent data representations consistently outperforms modelling on the original data. The results show an overall accuracy of 96.12%, whereas additional models achieve an average precision rate of 95.12% and an average recall rate of 94.40%. The semi-supervised learn- ing approaches in this thesis compare favourable to scANVI in terms of accuracy, average precision rate, average recall rate and average f1-score. Moreover, results for alternative scenarios, in which cell types among training and test data do not perfectly overlap, are reported in this thesis.

## Methodology framework
![alt text](https://github.com/Thijsq/thesis/blob/master/Methodology.png)

The methodology consists of seven different steps, with corresponding code files:

### Code:
#### Data preprocessing

#### Grouped LASSO in original data

#### Self-learning grouped LASSO in original data

#### scVI - transition to latent dimensionality

#### Grouped LASSO in latent space

#### Self-learning grouped LASSO in latent space

#### Label propagation in latent space

#### Label spreading in latent space

#### Classification reports
