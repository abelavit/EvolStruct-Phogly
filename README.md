# EvolStruct-Phogly

Matlab Codes Usage

The matlab codes provided in this file can be used for two purposes. The details of the usage is given below.

Option 1: Executing the Train and Test datasets used for building the PhoglyStruct predictor

For this option, please run the m files named EvolStruct_Phogly_option1, followed by CKSAAP and Phogly_PseAAC. In this scenario, EvolStruct_Phogly_option1 code carries out the 10 fold cross-validation on the LibSVM classifier using the provided ‘Test_Sets’ datasets which are the same datasets used for building the EvolStruct-Phogly predictor. CKSAAP and Phogly_PseAAC codes return the performance metrics against the same test dataset. The performance metrics for the methods EvolStruct-Phogly, CKSAAP_PhoglySite and Phogly_PseAAC can be viewed under the variables Results_Avg, Results_Avg_CKSAAP and Results_Avg_Phogly_PseAAC respectively after the code execution. Furthermore, the AUC of EvolStruct-Phogly and CKSAAP_PhoglySite methods can be accessed under variable names AUC and AUC2 respectively.

Option 2: Generating new Train and Test datasets for PhoglyStruct method

For this option, please run the m files named EvolStruct_Phogly_Option2, followed by CKSAAP and Phogly_PseAAC. In this scenario, EvolStruct_Phogly_Option2 code extracts the structural properties and evolutionary information of each lysine k and carries out the 10 fold cross-validation on the LibSVM classifier using the newly generated dataset which very likely would be different compared to datasets used for building the EvolStruct-Phogly predictor. The ‘Test_Sets’ in option 1 was obtained the same way as the test sets here. CKSAAP and Phogly_PseAAC codes return the performance metrics against the newly generated test dataset. The performance metrics for the methods EvolStruct-Phogly, CKSAAP_PhoglySite and Phogly_PseAAC can be viewed under the variables Results_Avg, Results_Avg_CKSAAP and Results_Avg_Phogly_PseAAC respectively after the code execution. Furthermore, the AUC of EvolStruct-Phogly and CKSAAP_PhoglySite methods can be accessed under variable names AUC and AUC2 respectively.

Note: Before executing the codes, please make sure to change the libsvm-weights-3.22 file directory in the codes to your own.

Footnotes:

To find in detail the CKSAAP_PhoglySite feature extraction method for each lysine k, please see the code file named ‘CKSAAP_Preprocessing’. After the code execution, features are saved in the Final_Data variable. This is the same file used when comparing for the CKSAAP_PhoglySite method in the above two options.

The file also contains FASTA format of the phosphoglycerylation dataset which was used to obtained the predictions of all lysine k from the Phogly–PseAAC webserver accessible at http://app.aporc.org/Phogly-PseAAC/
