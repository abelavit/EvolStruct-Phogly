
% Compare the samples in testing dataset with already obtained labels from
% Phogly–PseAAC webserver at http://app.aporc.org/Phogly-PseAAC/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load Phogly_PseAAC_Result % Contains the prediction by the Phogly–PseAAC predictor

Whole_Test = {'Fold1_10' 'Fold2_10' 'Fold3_10' 'Fold4_10' 'Fold5_10' 'Fold6_10' 'Fold7_10' 'Fold8_10' 'Fold9_10' 'Fold10_10'}; 

for k = 1:10

    st = num2str(k);
    test_name = strcat('Fold', st, '_10');
    fold_file = [eval(test_name)];
    
    FN = 0;
    FP = 0;
    TN = 0;
    TP = 0;

    for i=1:size(fold_file,1)

        phosphoglycerylation_match = 0;  

        for j=1:size(result_file,1)

           if (strcmp(fold_file{i,1}, result_file{j,2}) == true) && (result_file{j,3} == fold_file{i,4})

                phosphoglycerylation_match = 1; 

                break % As soon as ks from fold_file and result_file match, exit the for loop

            end

        end
        if (fold_file{i,3} == '1') && (phosphoglycerylation_match == 1)

            TP = TP + 1;

        elseif (fold_file{i,3} == '1') && (phosphoglycerylation_match == 0)

            FN = FN + 1;

        elseif (fold_file{i,3} == '0') && (phosphoglycerylation_match == 1)

            FP = FP + 1;

        elseif (fold_file{i,3} == '0') && (phosphoglycerylation_match == 0)

            TN = TN + 1;

        end

    end
    sen = TP/(TP+FN);
    spe = TN/(TN+FP);
    pre = TP/(TP+FP);
    accuracy = (TN+TP)/(FN+FP+TN+TP);
    mcc = ((TN*TP)-(FN*FP))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    Results_Phogly_PseAAC(1,k) = sen; 
    Results_Phogly_PseAAC(2,k) = spe; 
    Results_Phogly_PseAAC(3,k) = pre; 
    Results_Phogly_PseAAC(4,k) = accuracy; 
    Results_Phogly_PseAAC(5,k) = mcc;
end
Results_Avg_Phogly_PseAAC = sum(Results_Phogly_PseAAC,2)/10;

% Display the result
metrics = {'sensitivity'; 'specificity'; 'precision'; 'accuracy'; 'mcc'};
Ten_Fold_CV_PseAAC = [metrics, num2cell(Results_Avg_Phogly_PseAAC)] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
