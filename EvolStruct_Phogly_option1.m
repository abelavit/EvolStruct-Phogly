

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
% Contruction of training sets
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load Test_Sets

Whole_Data = {'Fold1_10' 'Fold2_10' 'Fold3_10' 'Fold4_10' 'Fold5_10' 'Fold6_10' 'Fold7_10' 'Fold8_10' 'Fold9_10' 'Fold10_10'}; 

for i = 1:10
     
    st = num2str(i);
    s = strcat('Fold', st, '_10');
    
    k = 1;
    for j = 1:10
        if regexp(Whole_Data{1,j}, s) == 1 % For each part, construct the corresponding unoptimized training set composed of other parts
            % Do nothing
        else
            index_train{1,k} = j;
            k = k+1;
        
        end
        
    end
    rt = num2str(index_train{1,1});
    r = strcat('Fold', rt, '_10');
    st = num2str(index_train{1,2});
    s = strcat('Fold', st, '_10');
    tt = num2str(index_train{1,3});
    t = strcat('Fold', tt, '_10');
    ut = num2str(index_train{1,4});
    u = strcat('Fold', ut, '_10');
    vt = num2str(index_train{1,5});
    v = strcat('Fold', vt, '_10');
    wt = num2str(index_train{1,6});
    w = strcat('Fold', wt, '_10');
    xt = num2str(index_train{1,7});
    x = strcat('Fold', xt, '_10');
    yt = num2str(index_train{1,8});
    y = strcat('Fold', yt, '_10');
    zt = num2str(index_train{1,9});
    z = strcat('Fold', zt, '_10');
    
    Train{i} = [eval(r); eval(s); eval(t); eval(u); eval(v); eval(w); eval(x); eval(y); eval(z)]; % Save the other 9 parts as unoptimized Training set
    
end

% Saving the training set into different variable

Train1 = Train{1};
Train2 = Train{2};
Train3 = Train{3};
Train4 = Train{4};
Train5 = Train{5};
Train6 = Train{6};
Train7 = Train{7};
Train8 = Train{8};
Train9 = Train{9};
Train10 = Train{10};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------
% Training the Classifier and obtaining performance metrics
%-----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd 'C:\Users\User\Desktop\libsvm-weights-3.22\matlab'

Whole_Test = {'Fold1_10' 'Fold2_10' 'Fold3_10' 'Fold4_10' 'Fold5_10' 'Fold6_10' 'Fold7_10' 'Fold8_10' 'Fold9_10' 'Fold10_10'}; 
Whole_Train = {'Train1' 'Train2' 'Train3' 'Train4' 'Train5' 'Train6' 'Train7' 'Train8' 'Train9' 'Train10'}; 

for i = 1:10

    clear Test_Data Train_Data test_label train_label Test_label Train_label Test Train
    st = num2str(i);
    test_name = strcat('Fold', st, '_10');
    train_name = strcat('Train', st);
    
    % Prepare the train, test, train label and test label sets
    Test_Data = [eval(test_name)];
    Train_Data = [eval(train_name)];
    
    Test = cell2mat(Test_Data(:,2));
    Train = cell2mat(Train_Data(:,2));

    test_label = cell2mat(Test_Data(:,3));
    train_label = cell2mat(Train_Data(:,3));
    Test_label = str2num(test_label);
    Train_label = str2num(train_label);

    % Train LibSVM-weights. Empty square brackets denote we are not
    % supplying weights
    model=svmtrain([],Train_label,Train,['-s 0 -t 1 -c 1 -g 1']);
    
    % Test the classifier
    [pred,acc,prob_values]=svmpredict(Test_label,Test,model);
    
    % Save the predicted labels and the true labels
    prediction{i} = pred;
    True_label{i} = Test_label;

    % Obtaining the FN, FP, TN and TP values
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
    
    % Calculating the performance metrics
    sen = TP/(TP+FN);
    spe = TN/(TN+FP);
    pre = TP/(TP+FP);
    accuracy = (TN+TP)/(FN+FP+TN+TP);
    mcc = ((TN*TP)-(FN*FP))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    Results(1,i) = sen; 
    Results(2,i) = spe; 
    Results(3,i) = pre; 
    Results(4,i) = accuracy; 
    Results(5,i) = mcc;
end
Results_Avg = sum(Results,2)/10; % 10 fold CV result 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----------------------------------------------------------------------
% AUC calculation and displaying of all the results
%-----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classes = [True_label{1}; True_label{2}; True_label{3}; True_label{4}; True_label{5}; True_label{6}; True_label{7}; True_label{8}; True_label{9}; True_label{10}];
scores = [prediction{1}; prediction{2}; prediction{3}; prediction{4}; prediction{5}; prediction{6}; prediction{7}; prediction{8}; prediction{9}; prediction{10}];

[X,Y,T,AUC] = perfcurve(classes,scores,1);

% Display the result
metrics = {'sensitivity'; 'specificity'; 'precision'; 'accuracy'; 'mcc'};
Ten_Fold_CV = [metrics, num2cell(Results_Avg)] 
AUC

