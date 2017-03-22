function [b_acc_opls] = Complete_Histograms( S,binsize,VariableNames,dataname,plotme)
%% Created by Thomas Vijverberg on 09-06-2015 at Radboud University Nijmegen
% Last edited by Thomas Vijverberg on 09-06-2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input data structure (STRUCT) 'S' with dimension 1 x (total_measurements)
    %Must contain atleast contain:
    %Data input (S.Data)
    %IDs input (S.ID)
    %Labels S.Labels
    
%Input binsize for histograms (INT) 20-1500 recommended

%Input VariableNames matrix (MAT) with dimension 1 x (number_variables)

%Input if you want the resulting DAMACY PLOT (BOOL) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preprocessing

%% Outlier removal 
% Calculate 95% quantiles (5% considered outlier) for each ID,variable.
% If cell of individual has 4 or more variables considered outlier, it is
% removed. This is done to correct instrumental errors in FACS data.

data = vertcat(S.Data)';
N_Var = length(S(1).Data(1,:));
N_ID = length(S);
paired = 0;

for i = 1:N_Var
    y(i,:) = quantile(data(i,:),[0.025 0.25 0.50 0.75 0.975]);
end
S_temp = S;
S_temp2 =S;
S_O = S;

for j = 1:N_ID
    for i = 1:N_Var
        S_temp(j).Data(:,i)=S(j).Data(:,i) <= y(i,1);
        S_temp2(j).Data(:,i)=S(j).Data(:,i) >= y(i,5);
    end
    S_temp(j).Data2(:,1) = ~((sum(S_temp(j).Data,2)+sum(S_temp2(j).Data,2)) > 2);
    S_O(j).Data = S(j).Data(S_temp(j).Data2,:);
end

%% Meancentering and scaling 
%MEANCENTERING
%1 mean center over all controls
%2 median center over all controls
%3 mean center over all individuals 
%4 median center over all individuals
% SCALING
%1 Pareto scaling complete dataset
%2 UV scaling over control group
%3 UV scaling over individuals


N_vars = length(S(1).Data(1,:));
N_id = length(S);
S_MC = unpaired_centering(S_O,3);
S_A = scaling(S_MC, vertcat(S.Labels), 1);

data = vertcat(S_A.Data);
minimum_value = min(data);
maximum_value = max(data);
%% Base Model
% Create histogram grid (edges) for each variable over complete dataset
% Length of grid will be binsize -2 (one missing due to matlab histcounts
% implementation, one missing for removing zero measurements in first bin)

for i = 1:length(minimum_value)
    edges(i,:) = linspace(minimum_value(i),maximum_value(i),binsize);
end

% Allocate memory for histogram structure
S_hist = struct('Data',[],'Labels',[],'ID',[]);
for i = 1:length(S)
    S_hist(i).Data = zeros(binsize-1,N_vars);
end
X = zeros(N_id,N_vars*(binsize-2));

% Create histograms from original data using matlab function HISTCOUNTS
% Smooth histograms using smoothn MATLAB FILE EXCHANGE:
for i = 1:length(S)
    for j = 1:length(S(1).Data(1,:))
        S_hist(i).Data(:,j) = histcounts(S_A(i).Data(:,j),edges(j,:));
        S_hist(i).Data(:,j) = smoothn(S_hist(i).Data(:,j),0.45,'robust');
    end
    S_hist(i).Data = S_hist(i).Data(2:binsize-1,:);
    X(i,:) = reshape(S_hist(i).Data,1,[]);
end

%% Top Model OPLS
% Preparation of data for OPLS function. Labels have to be 1 and -1

ID = vertcat(S.ID);
Labels = vertcat(S.Labels);
Y2 = Labels == 0;
Y1 = ~Y2;
Y = Y1 - Y2;
mean_Y = mean(Y);
X = X - repmat(mean(X), size(X,1), 1);
Ym = Y - mean(Y);

% Perform OPLS regression and classification
%b_acc is the OPLS accuracy calculated with either LOOCV / LTOCV
%w_opls is the OPLS weights vector, the regression model
%yhat are the predicted labels using LOOCV / LTOCV classification
[b_acc_opls,n_LV,w_opls,yhat_opls] = DAMACY_top(X,ID,Y,5,paired);

%% Plot
if plotme == 1
    Analysis_Result_1D(S_A,w_opls,edges,binsize,b_acc_opls,yhat_opls,Y,VariableNames,'OPLS',dataname,paired);
end

%% Top Model SVM
% SVM_model = svmtrain(Labels,X,'-c 700 -t 2 -g 2 -r 1 -q');
% [yhatt_svm, b_acc_svm, yhat_svm] = svmpredict(Labels, X, SVM_model, '');
% [b_acc] = DAMACY_SVM_top(X,ID,Labels);
% w_svm = SVM_model.SVs' * SVM_model.sv_coef;
% Analysis_Result_1D(w_svm,edges,binsize,b_acc_svm(1)/100,yhat_svm,VariableNames, 'SVM',dataname);
end

