
% Compare and extract similar Testing data from CKSAAP_Data. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load CKSAAP_Data % This is the dataset containing CKSAAP features for each lysine K
load Feature_Label
clear Train

Whole_Data = {'Fold1_10' 'Fold2_10' 'Fold3_10' 'Fold4_10' 'Fold5_10' 'Fold6_10' 'Fold7_10' 'Fold8_10' 'Fold9_10' 'Fold10_10'}; 

%Test Sets
for h = 1:10
     
    st = num2str(h);
    s = strcat('Fold', st, '_10');
    fold_file = [eval(s)];
    
    for i=1:size(fold_file,1)
    
    for j=1:size(Final_Data,1)

        if (strcmp(Final_Data{j,1}, fold_file{i,1}) == true) && (Final_Data{j,4} == fold_file{i,4})

                   NewFold(i,:) = Final_Data(j,:);
                   
                   break % As soon as ks from fold_file and result_file match, exit the for loop
        end
    end
 
    end
    
    Fold{h} = NewFold; % Save each generated fold under one variable
    
end

%Train Sets
for k = 1:10
    z = 1;
    for l = 1:10
    
        if k == l
            % Do nothing
        else
            
            T(1,z) = l;
            z = z + 1;
        end
        
    end

    Train{k} = [Fold{T(1,1)}; Fold{T(1,2)}; Fold{T(1,3)}; Fold{T(1,4)}; Fold{T(1,5)}; Fold{T(1,6)}; Fold{T(1,7)}; Fold{T(1,8)}; Fold{T(1,9)}]; % use string name to save cell array of same name into Data
    
end

%F_Score calculation for each training set
for m = 1:10
    % Number of +ve and -ve samples and Separating +ve and -ve samples
    % The if conditions below will result in error if sample label is of type
    % double
    Data2 = Train{m};

    Num_of_Pos  = 0;
    Num_of_Neg  = 0;
    for i = 1:size(Data2,1)
        if Data2{i,3} == '1'
           Num_of_Pos = Num_of_Pos + 1;
           Pos_sample(Num_of_Pos,:) = Data2(i,:); 
        end
        if Data2{i,3} == '0'
           Num_of_Neg = Num_of_Neg + 1;
           Neg_sample(Num_of_Neg,:) = Data2(i,:);
        end    
    end
    % Sum the columns of the Final_Data
    for i = 5:size(Data2,2)

        Summ = 0;
        for j = 1:size(Data2,1)

          Summ = Summ + Data2{j,i};  

        end
        Feature_sum{1,i} = Summ;
    end

    % Mean j-th feature in whole
    for i = 5:size(Feature_sum, 2)

         Feature_mean_whole{1,i} = Feature_sum{1,i}/size(Data2,1);    

    end

    % Mean j-th feature in positive

    % Sum the columns of the Pos_sample
    for i = 5:size(Pos_sample,2)

        Summ = 0;
        for j = 1:size(Pos_sample,1)

          Summ = Summ + Pos_sample{j,i};  

        end
        Feature_sum_pos{1,i} = Summ;
    end
    % Mean j-th feature in Pos_sample
    for i = 5:size(Feature_sum_pos, 2)

         Feature_mean_pos{1,i} = Feature_sum_pos{1,i}/size(Pos_sample,1);    

    end

    % Mean j-th feature in negative

    % Sum the columns of the Neg_sample
    for i = 5:size(Neg_sample,2)

        Summ = 0;
        for j = 1:size(Neg_sample,1)

          Summ = Summ + Neg_sample{j,i};  

        end
        Feature_sum_neg{1,i} = Summ;
    end
    % Mean j-th feature in Neg_sample
    for i = 5:size(Feature_sum_neg, 2)

         Feature_mean_neg{1,i} = Feature_sum_neg{1,i}/size(Neg_sample,1);    

    end

    % Denominator part of F_Score formula

    % Pos part of deno
    for i = 5:size(Pos_sample, 2)
        Summ = 0;
        for j = 1:size(Pos_sample, 1)

            Summ = Summ + (Pos_sample{j,i} - Feature_mean_pos{1,i})^2;

        end
        Pos_deno{1,i} = (1/(size(Pos_sample, 1) - 1)) * Summ;
    end
    % Neg part of deno
    for i = 5:size(Neg_sample, 2)
        Summ = 0;
        for j = 1:size(Neg_sample, 1)

            Summ = Summ + (Neg_sample{j,i} - Feature_mean_neg{1,i})^2;

        end
        Neg_deno{1,i} = (1/(size(Neg_sample, 1) - 1)) * Summ;
    end

    % F-Score Calculation

    for i = 5:size(Data2, 2)

        if (Pos_deno{1,i} + Neg_deno{1,i}) == 0
            F_Score{1,i} = 0;

        else
            F_Score{1,i} = ((Feature_mean_pos{1,i}-Feature_mean_whole{1,i})^2 + (Feature_mean_neg{1,i}-Feature_mean_whole{1,i})^2)/(Pos_deno{1,i} + Neg_deno{1,i});
        end
    end

    Fscore_n_FeatureLabel = [F_Score; feature_label];

    A = Fscore_n_FeatureLabel';
    B = A(5:size(A,1),:);
    C = sortrows(B,1);
    Ju_Score{m} = flipud(C);
    
    clear Pos_sample;
    clear Neg_sample;
    clear Feature_sum;
    clear Feature_mean_whole;
    clear Feature_sum_pos;
    clear Feature_mean_pos;
    clear Feature_sum_neg;
    clear Feature_mean_neg;
    clear Pos_deno;
    clear Neg_deno;
    clear F_Score;
    clear A;
    clear B;
    clear C;
    
end

% Select first 300 features of both train and test sets using F_score (Ju's score)

for n = 1:10
    
    Data = Train{n};
    Data_300 = Data(:,1:4);

    for i=1:300

        for j=5:size(feature_label, 2)

            if strcmp(feature_label(1,j),Ju_Score{n}(i,2)) == true

                Data_300(:,i+4) = Data(:,j);

            end
        end

    end
    Final_Train{n} = Data_300;    
    
    clear Data_300; 
    clear Data;
    
    Data = Fold{n};
    Data_300 = Data(:,1:4);

    for i=1:300

        for j=5:size(feature_label, 2)

            if strcmp(feature_label(1,j),Ju_Score{n}(i,2)) == true

                Data_300(:,i+4) = Data(:,j);

            end
        end

    end
    Final_Fold{n} = Data_300;
    
    clear Data_300; 
    clear Data;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Training the Classifier and obtaining performance metrics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd 'C:\Users\User\Desktop\libsvm-weights-3.22\matlab'

for i = 1:10
    
    clear Test_Data Train_Data test_label train_label Test_label Train_label Test Train
    % Prepare the train, test, train label and test label sets
    Test_Data = Final_Fold{i};
    Train_Data = Final_Train{i};

    for k = 1:size(Test_Data,1)
        Test(k,:) = cell2mat(Test_Data(k,5:size(Test_Data,2)));
    end
    for k = 1:size(Train_Data,1)
        Train(k,:) = cell2mat(Train_Data(k,5:size(Train_Data,2)));
    end

    test_label = cell2mat(Test_Data(:,3));
    train_label = cell2mat(Train_Data(:,3));
    Test_label = str2num(test_label);
    Train_label = str2num(train_label);
    
    model=svmtrain([],Train_label,Train,['-s 0 -t 1 -c 1 -g 1']);

    [pred,acc,prob_values]=svmpredict(Test_label,Test,model);
    
    predicted_Labels{i} = pred;
    True_labels{i} = Test_label;

    FN = 0;
    FP = 0;
    TN = 0;
    TP = 0;
    for j = 1:size(Test_label,1)
        if Test_label(j) == 1
            if pred(j) == 1
                TP = TP + 1;
            else
                FN = FN + 1; 
            end
        else
            if pred(j) == 1
                FP = FP + 1;
            else
                TN = TN + 1; 
            end
        end
    end
    
    if (TP+FN) == 0
        sen = 0;
    else
        sen = TP/(TP+FN);
    end
    if (TN+FP) == 0
        spe = 0;
    else
        spe = TN/(TN+FP);
    end
    if (TP+FP) == 0
        pre = 0;
    else
        pre = TP/(TP+FP);
    end     
    accuracy = (TN+TP)/(FN+FP+TN+TP);
    if ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) == 0
        mcc = 0;
    else
        mcc = ((TN*TP)-(FN*FP))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end    
    
    Results_CKSAAP(1,i) = sen; 
    Results_CKSAAP(2,i) = spe; 
    Results_CKSAAP(3,i) = pre; 
    Results_CKSAAP(4,i) = accuracy; 
    Results_CKSAAP(5,i) = mcc;
end
Results_Avg_CKSAAP = sum(Results_CKSAAP,2)/10;

% AUC calculation
classes = [True_labels{1}; True_labels{2}; True_labels{3}; True_labels{4}; True_labels{5}; True_labels{6}; True_labels{7}; True_labels{8}; True_labels{9}; True_labels{10}];
scores = [predicted_Labels{1}; predicted_Labels{2}; predicted_Labels{3}; predicted_Labels{4}; predicted_Labels{5}; predicted_Labels{6}; predicted_Labels{7}; predicted_Labels{8}; predicted_Labels{9}; predicted_Labels{10}];

[X2,Y2,T2,AUC2] = perfcurve(classes,scores,1);

% Display the result
metrics = {'sensitivity'; 'specificity'; 'precision'; 'accuracy'; 'mcc'};
Ten_Fold_CV_CKSAAP = [metrics, num2cell(Results_Avg_CKSAAP)] 
AUC2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
