# predict_HPC_AEM

1. Files with df_predict_pc_composite_dopant#0.csv are the predicted HPC values of oxides containing 8 dopants with varying doping fraction, and the doping fraction are 10%, 20%, 30%, 40%, and 50%.
2. File named as "Original dataset for training models_with_References.xls" are the original data containing 784 perovskite oxides and their corresponding source references.
3. Script named as "XGB_github.R " provides the details to train machine learning model with XGBoost algorithm.
4. Script named as "RF_github.R " provides the details to train machine learning model with random forest algorithm.
5. Script named as "Asite_predict_CoFe.R" provides the details to predict HPC of AB_{1-x}B1_{x}O3-type oxides.
6. Script named as "AA1B_descriptor_customize.R" provides the details to customize the descriptors, which would be used for the following model training and prdiction.
7. File named as "elements_data.csv" provides all the quantative attributes to represent alphabetic element name.
8. File named as "saveXGBmodel_rep10cv10_20231024.rdata" is an intermediate process file, which are saved from model training and would be used for model predicting.
9. The zip file named as "Examples of DFT calculation files.zip" contains all the source files for DFT calculation.
   
