
%-------------------------------------------------------------
% First Phase: Data Extraction and construction of testing set
%-------------------------------------------------------------
load Phosphoglycerylationstruct % Load the phosphoglycerylation dataset (raw file from which we will extract the data from)

% Data extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Unprocessed_data = DB_Phosphoglycerylation;

Field = size(Unprocessed_data,2); % Columns of unprocessed data. It is number of protein sequences 

z=0;
for a = 1:Field
    z = z + size(strfind(Unprocessed_data(a).seq{1},'K'), 2); % Finding total lysine K in the entire protein sequences
end

Final_Data = cell(z,4); % Create data file which will contain info of each lysine k
k=0;  

for l = 1:Field
    K_locations_field{l} = strfind(Unprocessed_data(l).seq{1},'K'); % Saves the locations of K found in protein sequences
end

window = 20; % Window of 20 indicates 20 upstream and 20 downstream of lysine k for Evolutionary information
window2 = 3; % Window2 of 3 indicates 3 upstream and 3 downstream of lysine k for Structural Properties

for m = 1:Field
    for n=1:size(K_locations_field{m},2) % loop in all the K locations
        
        k = k+1; % Increment to next data saving location
        Final_Data{k,1} = Unprocessed_data(m).name; % Save protein name at 1st position
        
        % Evolutionary information (PSSM)
        if K_locations_field{m}(n) <= window % Location K close to N terminus
            matrix_a = Unprocessed_data(m).pssm_prob(1:K_locations_field{m}(n)+ window,:); % Matrix containing pssm_prob from start of protein seq till 20 downstream of K
            matrix_b = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+ window,:); % Matrix containing pssm_prob from location 2x the location of K till end of downstream 
            matrix_c = flipud(matrix_b); % Mirroring the matrix_b (flipping on horizontal axis) 
            Final_Data{k,5} = [matrix_c; matrix_a]; % Concatenate the two matrices which contains sufficient upstream and downstream amino acid. This is the matrix M 
       
        elseif K_locations_field{m}(n) > (Unprocessed_data(m).len - window) % Location K close to C terminus
            matrix_d = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)-window:Unprocessed_data(m).len,:); % Matrix containing pssm_prob from 20 downstream of K till end of protein seq
            matrix_e = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)-window:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Matrix containing pssm_prob upstream of K till the value to be mirrored
            matrix_f = flipud(matrix_e); % Mirroring the matrix_e (flipping on horizontal axis)
            Final_Data{k,5} = [matrix_d; matrix_f]; % Concatenate the two matrices which contains sufficient upstream and downstream amino acid. This is the matrix M
            
        else % Location is good and does not require mirroring
            Final_Data{k,5} = Unprocessed_data(m).pssm_prob(K_locations_field{m}(n)-window:K_locations_field{m}(n)+window,:); % Save pssm_prob matrix (20 up and 20 down stream of K). This is the matrix M
        
        end
        
        % Structural properties (SPpre)
        if K_locations_field{m}(n) <= window2 % Location K close to N terminus

            matrix_a_1 = Unprocessed_data(m).ASA(1:K_locations_field{m}(n)+window2,:); % Vector containing ASA from start of protein seq till 3 downstream of K
            matrix_b_1 = Unprocessed_data(m).ASA(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+window2,:); % % Vector containing ASA from location 2x the location of K till end of downstream 
            matrix_c_1 = flipud(matrix_b_1); % Mirroring the matrix_b_1 (flipping on horizontal axis) 
            Final_1_1 = [matrix_c_1; matrix_a_1]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid 

            matrix_a_2 = Unprocessed_data(m).SSpre(1:K_locations_field{m}(n)+window2,:); % Matrix containing SSpre from start of protein seq till 3 downstream of K. SSpre contains prediction of amino acid to one of the three conformations helix, strand and coil
            matrix_b_2 = Unprocessed_data(m).SSpre(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+window2,:); % Matrix containing SSpre from location 2x the location of K till end of downstream 
            matrix_c_2 = flipud(matrix_b_2); % Mirroring the matrix_b_2 (flipping on horizontal axis)
            Final_1_2 = [matrix_c_2; matrix_a_2]; % Concatenate the two matrices which contains sufficient upstream and downstream amino acid
            
            matrix_a_3 = Unprocessed_data(m).Phi(1:K_locations_field{m}(n)+window2,:); % Vector containing Phi from start of protein seq till 3 downstream of K
            matrix_b_3 = Unprocessed_data(m).Phi(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+window2,:); % Vector containing Phi from location 2x the location of K till end of downstream 
            matrix_c_3 = flipud(matrix_b_3); % Mirroring the matrix_b_3 (flipping on horizontal axis) 
            Final_1_3 = [matrix_c_3; matrix_a_3]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid 

            matrix_a_4 = Unprocessed_data(m).Psi(1:K_locations_field{m}(n)+window2,:); % Vector containing Psi from start of protein seq till 3 downstream of K
            matrix_b_4 = Unprocessed_data(m).Psi(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+window2,:); % Vector containing Psi from location 2x the location of K till end of downstream 
            matrix_c_4 = flipud(matrix_b_4); % Mirroring the matrix_b_4 (flipping on horizontal axis)  
            Final_1_4 = [matrix_c_4; matrix_a_4]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            matrix_a_5 = Unprocessed_data(m).Theta(1:K_locations_field{m}(n)+window2,:); % Vector containing Theta from start of protein seq till 3 downstream of K
            matrix_b_5 = Unprocessed_data(m).Theta(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+window2,:); % Vector containing Theta from location 2x the location of K till end of downstream 
            matrix_c_5 = flipud(matrix_b_5); % Mirroring the matrix_b_5 (flipping on horizontal axis) 
            Final_1_5 = [matrix_c_5; matrix_a_5]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid 

            matrix_a_6 = Unprocessed_data(m).Tau(1:K_locations_field{m}(n)+window2,:); % Vector containing Tau from start of protein seq till 3 downstream of K
            matrix_b_6 = Unprocessed_data(m).Tau(K_locations_field{m}(n)+K_locations_field{m}(n):K_locations_field{m}(n)+window2,:); % Vector containing Tau from location 2x the location of K till end of downstream 
            matrix_c_6 = flipud(matrix_b_6); % Mirroring the matrix_b_6 (flipping on horizontal axis) 
            Final_1_6 = [matrix_c_6; matrix_a_6]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            F1 = transpose([Final_1_1 Final_1_2 Final_1_3 Final_1_4 Final_1_5 Final_1_6]); % Concatenate the matrices containing the 8 properties of each amino acid
            Final_Data{k,6} = F1(:)'; % Transpose the concatenated matrix and save it. This makes up the features for Lysine K close to the start of the protein sequence. This is SPpre 
            
        elseif K_locations_field{m}(n) > (Unprocessed_data(m).len - window2) % Location K close to C terminus

            matrix_d_1 = Unprocessed_data(m).ASA(K_locations_field{m}(n)-window2:Unprocessed_data(m).len,:); % Vector containing ASA from 3 upstream of K till end of protein seq
            matrix_e_1 = Unprocessed_data(m).ASA(K_locations_field{m}(n)-window2:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Vector containing ASA from 3 upstream of K till the value to be mirrored
            matrix_f_1 = flipud(matrix_e_1); % Mirroring the matrix_e_1 (flipping on horizontal axis)
            Final_2_1 = [matrix_d_1; matrix_f_1]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            matrix_d_2 = Unprocessed_data(m).SSpre(K_locations_field{m}(n)-window2:Unprocessed_data(m).len,:); % Matrix containing SSpre from 3 upstream of K till end of protein seq
            matrix_e_2 = Unprocessed_data(m).SSpre(K_locations_field{m}(n)-window2:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Matrix containing SSpre from 3 upstream of K till the value to be mirrored
            matrix_f_2 = flipud(matrix_e_2); % Mirroring the matrix_e_2 (flipping on horizontal axis)
            Final_2_2 = [matrix_d_2; matrix_f_2]; % Concatenate the two matrices which contains sufficient upstream and downstream amino acid

            matrix_d_3 = Unprocessed_data(m).Phi(K_locations_field{m}(n)-window2:Unprocessed_data(m).len,:); % Vector containing Phi from 3 upstream of K till end of protein seq
            matrix_e_3 = Unprocessed_data(m).Phi(K_locations_field{m}(n)-window2:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Vector containing Phi from 3 upstream of K till the value to be mirrored
            matrix_f_3 = flipud(matrix_e_3); % Mirroring the matrix_e_3 (flipping on horizontal axis)
            Final_2_3 = [matrix_d_3; matrix_f_3]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            matrix_d_4 = Unprocessed_data(m).Psi(K_locations_field{m}(n)-window2:Unprocessed_data(m).len,:); % Vector containing Psi from 3 upstream of K till end of protein seq
            matrix_e_4 = Unprocessed_data(m).Psi(K_locations_field{m}(n)-window2:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Vector containing Psi from 3 upstream of K till the value to be mirrored
            matrix_f_4 = flipud(matrix_e_4); % Mirroring the matrix_e_4 (flipping on horizontal axis)
            Final_2_4 = [matrix_d_4; matrix_f_4]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            matrix_d_5 = Unprocessed_data(m).Theta(K_locations_field{m}(n)-window2:Unprocessed_data(m).len,:); % Vector containing Theta from 3 upstream of K till end of protein seq
            matrix_e_5 = Unprocessed_data(m).Theta(K_locations_field{m}(n)-window2:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Vector containing Theta from 3 upstream of K till the value to be mirrored
            matrix_f_5 = flipud(matrix_e_5); % Mirroring the matrix_e_5 (flipping on horizontal axis)
            Final_2_5 = [matrix_d_5; matrix_f_5]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            matrix_d_6 = Unprocessed_data(m).Tau(K_locations_field{m}(n)-window2:Unprocessed_data(m).len,:); % Vector containing Tau from 3 upstream of K till end of protein seq
            matrix_e_6 = Unprocessed_data(m).Tau(K_locations_field{m}(n)-window2:Unprocessed_data(m).len-(((Unprocessed_data(m).len-K_locations_field{m}(n))*2)+1),:); % Vector containing Tau from 3 upstream of K till the value to be mirrored
            matrix_f_6 = flipud(matrix_e_6); % Mirroring the matrix_e_6 (flipping on horizontal axis)
            Final_2_6 = [matrix_d_6; matrix_f_6]; % Concatenate the two vectors which contains sufficient upstream and downstream amino acid

            F1 = transpose([Final_2_1 Final_2_2 Final_2_3 Final_2_4 Final_2_5 Final_2_6]); % Concatenate the matrices containing the 8 properties of each amino acid
            Final_Data{k,6} = F1(:)'; % Transpose the concatenated matrix and save it. This makes up the features for Lysine K close to the end of the protein sequence. This is SPpre 
            
        else % Location is good and does not require mirroring
            Final_3_1 = Unprocessed_data(m).ASA(K_locations_field{m}(n)-window2:K_locations_field{m}(n)+window2,:); % ASA vector (3 upstream and 3 downstream of K)
            Final_3_2 = Unprocessed_data(m).SSpre(K_locations_field{m}(n)-window2:K_locations_field{m}(n)+window2,:); % SSpre matrix (3 upstream and 3 downstream of K)
            Final_3_3 = Unprocessed_data(m).Phi(K_locations_field{m}(n)-window2:K_locations_field{m}(n)+window2,:); % Phi vector (3 upstream and 3 downstream of K)
            Final_3_4 = Unprocessed_data(m).Psi(K_locations_field{m}(n)-window2:K_locations_field{m}(n)+window2,:); % Psi vector (3 upstream and 3 downstream of K)
            Final_3_5 = Unprocessed_data(m).Theta(K_locations_field{m}(n)-window2:K_locations_field{m}(n)+window2,:); % Theta vector (3 upstream and 3 downstream of K)
            Final_3_6 = Unprocessed_data(m).Tau(K_locations_field{m}(n)-window2:K_locations_field{m}(n)+window2,:); % Tau vector (3 upstream and 3 downstream of K)

            F1 = transpose([Final_3_1 Final_3_2 Final_3_3 Final_3_4 Final_3_5 Final_3_6]); % Concatenate the matrices containing the 8 properties of each amino acid
            Final_Data{k,6} = F1(:)'; % Transpose the concatenated matrix and save it. This makes up the features for Lysine K which did not require mirroring. This is SPpre 
            
        end
        
        Final_Data{k,3} = Unprocessed_data(m).label{1}(K_locations_field{m}(n)); % Save class label
        Final_Data{k,4} = K_locations_field{m}(n); % Save the K's location in the protein sequence
    
    end
end

Training_Data1 = Final_Data; 

% PSSM + Bigram

for l=1:size(Training_Data1,1)
    Training_Data1{l,5} = Training_Data1{l,5}/100; % Divide pssm_prob by 100
    
    Bigram_Mat = zeros(size(Training_Data1{l,5}, 2),size(Training_Data1{l,5}, 2));
    for M=1:size(Training_Data1{l,5}, 2)
       for N=1:size(Training_Data1{l,5}, 2)
           for i=1:size(Training_Data1{l,5},1)-1
            Bigram_Mat(M,N) = Bigram_Mat(M,N) + Training_Data1{l,5}(i,M)*Training_Data1{l,5}(i+1,N); % Bigram calculation
           end
       end
    end
    transposed_Bigram_Mat = Bigram_Mat'; % Transpose Bigram_Mat
    Training_Data1{l,2} = [transposed_Bigram_Mat(:)' Final_Data{l,6}]; % Combine the PSSM+bigram and SPpre. This gives Feature Vector
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Separation of neg and pos samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Num_of_Pos  = 0;
Num_of_Neg  = 0;

for i = 1:size(Training_Data1,1)
    
    if Training_Data1{i,3} == '1'
       Num_of_Pos = Num_of_Pos + 1;
       Final_Pos(Num_of_Pos,:) = Training_Data1(i,:); 
    end
    if Training_Data1{i,3} == '0'
       Num_of_Neg = Num_of_Neg + 1;
       Neg_Samples(Num_of_Neg,:) = Training_Data1(i,:);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filtering of neg samples to bring class imbalance from 1:29 to 1:3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Filtered_Neg = Neg_Samples;

for a=1:size(Neg_Samples,1)
    
    for b=1:size(Training_Data1, 1)

        Distance = pdist2(Neg_Samples{a,2}, Training_Data1{b,2}); % Euclidean distance between the neg sample and a sample in whole data
   
        D(b, 1) = Distance; % Save the 3360 distances of entire data samples into 3360x1 matrix
        
    end
    
    for c=1:80 % Number of nearest neighbors (including itself as it also calculated from itself at one point in time) to be taken into account for filtering. Number of neighbors is actually 79 
   
        [M,I] = min(D); % Find minimum from the 3360x1 matix and also its index
        
        if Training_Data1{I,3} == '1' % If the minimum distance is for a +ve sample then delete the neg sample
            
            Filtered_Neg{a,2} = [];
            break
        end
       
        D(I) = [max(D)]; % Replace the minimum value with the maximum value of the matrix so that it wont point the same index as minimum in the next loop
        c = c+1;
    end
end

% Construct the neg samples after filtering
e = 1;
for j=1:size(Filtered_Neg,1)
    if isempty(Filtered_Neg{j,2}) == false
  
        Final_Neg(e,:) = Filtered_Neg(j,:); 
  
        e = e+1;
        
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dividing data into 10 parts for 10 fold cross validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Filtered_Data_SamplesSep = [Final_Pos; Final_Neg]; % Combine the final sets of positive and negatives samples after the filtering process

len_of_data = size(Filtered_Data_SamplesSep, 1);

extra = mod(len_of_data,10);
base_size = (len_of_data - extra)/10;

for i = 1:10
    if (i>(10-extra))
        fold_size(i) = base_size + 1;
        
    else
        fold_size(i) = base_size;
        
    end  
end

default_gap = [5 10 15 20 25 30 35 40 45 50]; % Randomly assigned numbers (Number of +ve samples in each parts)
loop = 1;

% Dividing the Dataset into 10 parts and making sure +ve and -ve samples
% are well distributed
while(loop == 1)

        Data_Shuffled1 = Filtered_Data_SamplesSep(randperm(len_of_data),:);
        Data_Shuffled = Data_Shuffled1(randperm(len_of_data),:);

        Fold1_10 = Data_Shuffled(1:fold_size(1),:); % First Part
        j1=0;
        for i=1:fold_size(1)
        if strfind(Fold1_10{i,3},'1') == true
        j1 = j1+1;
        end
        end
        default_gap(1) = j1;

        Fold2_10 = Data_Shuffled(fold_size(1)+1:sum(fold_size(1:2)),:); % Second Part
        j2=0;
        for i=1:fold_size(2)
        if strfind(Fold2_10{i,3},'1') == true
        j2 = j2+1;
        end
        end
        default_gap(2) = j2;

        Fold3_10 = Data_Shuffled(sum(fold_size(1:2))+1:sum(fold_size(1:3)),:); % Third Part
        j3=0;
        for i=1:fold_size(3)
        if strfind(Fold3_10{i,3},'1') == true
        j3 = j3+1;
        end
        end
        default_gap(3) = j3;

        Fold4_10 = Data_Shuffled(sum(fold_size(1:3))+1:sum(fold_size(1:4)),:); % Fourth Part
        j4=0;
        for i=1:fold_size(4)
        if strfind(Fold4_10{i,3},'1') == true
        j4 = j4+1;
        end
        end
        default_gap(4) = j4;

        Fold5_10 = Data_Shuffled(sum(fold_size(1:4))+1:sum(fold_size(1:5)),:); % Fifth Part
        j5=0;
        for i=1:fold_size(5)
        if strfind(Fold5_10{i,3},'1') == true
        j5 = j5+1;
        end
        end
        default_gap(5) = j5;

        Fold6_10 = Data_Shuffled(sum(fold_size(1:5))+1:sum(fold_size(1:6)),:); % Sixth Part
        j6=0;
        for i=1:fold_size(6)
        if strfind(Fold6_10{i,3},'1') == true
        j6 = j6+1;
        end
        end
        default_gap(6) = j6;

        Fold7_10 = Data_Shuffled(sum(fold_size(1:6))+1:sum(fold_size(1:7)),:); % Seventh Part
        j7=0;
        for i=1:fold_size(7)
        if strfind(Fold7_10{i,3},'1') == true
        j7 = j7+1;
        end
        end
        default_gap(7) = j7;

        Fold8_10 = Data_Shuffled(sum(fold_size(1:7))+1:sum(fold_size(1:8)),:); % Eighth Part
        j8=0;
        for i=1:fold_size(8)
        if strfind(Fold8_10{i,3},'1') == true
        j8 = j8+1;
        end
        end
        default_gap(8) = j8;

        Fold9_10 = Data_Shuffled(sum(fold_size(1:8))+1:sum(fold_size(1:9)),:); % Ninth Part
        j9=0;
        for i=1:fold_size(9)
        if strfind(Fold9_10{i,3},'1') == true
        j9 = j9+1;
        end
        end
        default_gap(9) = j9;

        Fold10_10 = Data_Shuffled(sum(fold_size(1:9))+1:sum(fold_size(1:10)),:); % Tenth Part
        j10=0;
        for i=1:fold_size(10)
        if strfind(Fold10_10{i,3},'1') == true
        j10 = j10+1;
        end
        end
        default_gap(10) = j10;
        
        smallest = min(default_gap); % Minimum number of +ve samples in the overall parts
        largest = max(default_gap); % Maximum number of +ve samples in the overall parts
        
        if ((largest - smallest) < 3) % If difference between minimum and maximum number of +ve samples is less then 3 then there is good distribution of +ves and -ves
            
            loop = 0;
            
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
% Second Phase: Contruction of training sets
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Third Phase: Training the Classifier and obtaining performance metrics
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
% Fourth Phase: AUC calculation and displaying of all the results
%-----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classes = [True_label{1}; True_label{2}; True_label{3}; True_label{4}; True_label{5}; True_label{6}; True_label{7}; True_label{8}; True_label{9}; True_label{10}];
scores = [prediction{1}; prediction{2}; prediction{3}; prediction{4}; prediction{5}; prediction{6}; prediction{7}; prediction{8}; prediction{9}; prediction{10}];

[X,Y,T,AUC] = perfcurve(classes,scores,1);

% Display the result
metrics = {'sensitivity'; 'specificity'; 'precision'; 'accuracy'; 'mcc'};
Ten_Fold_CV = [metrics, num2cell(Results_Avg)] 
AUC
