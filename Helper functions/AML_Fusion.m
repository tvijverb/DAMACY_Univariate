function [ b_acc_opls ] = AML_Fusion( S,binsize,VariableNames,dataname,plotme)
%AML_FUSION Summary of this function goes here
%   Detailed explanation goes here
ST= S;

Stemp = Get_Tube(ST,2);
N_ID = length(Stemp);
N_Var = length(Stemp(1).Data(1,:));
X = zeros(N_ID,N_Var*(binsize-2));
X2 = zeros(N_ID,N_Var*(binsize-2));
clear Stemp

for z = 1:2
    if z == 1
        S = Get_Tube(ST,2);
    else
        S = Get_Tube(ST,6);
    end
%% Preprocessing
%outlier removal
data = vertcat(S.Data)';

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
    S_temp(j).Data2(:,1) = ~((sum(S_temp(j).Data,2)+sum(S_temp2(j).Data,2)) > 3);
    S_O(j).Data = S(j).Data(S_temp(j).Data2,:);
end

% Scaling
N_vars = length(S(1).Data(1,:));
N_id = length(S);
S_MC = unpaired_centering(S,3);
S_A = scaling(S_MC, vertcat(S.Labels), 1);

data = vertcat(S_A.Data);
minimum_value = min(data);
maximum_value = max(data);
%% Base Model
for i = 1:length(minimum_value)
    edges(i,:) = linspace(minimum_value(i),maximum_value(i),binsize);
end

S_hist = struct('Data',[],'Labels',[],'ID',[]);
for i = 1:length(S)
    S_hist(i).Data = zeros(binsize-1,N_vars);
end

for i = 1:length(S)
    for j = 1:length(S(1).Data(1,:))
        S_hist(i).Data(:,j) = histcounts(S_A(i).Data(:,j),edges(j,:));
        S_hist(i).Data(:,j) = smoothn(S_hist(i).Data(:,j),0.45,'robust');
    end
    S_hist(i).Data = S_hist(i).Data(2:binsize-1,:);
    if z == 1
        X(i,:) = reshape(S_hist(i).Data,1,[]);
    else
        X2(i,:)= reshape(S_hist(i).Data,1,[]);
    end
end
if z == 1
    edgesz = edges;
else
    edges = vertcat(edgesz,edges);
end
end
X = horzcat(X,X2);
% ind = std(X) > 10^-6;
% X = X(:,ind);

%% Top Model OPLS
ID = vertcat(S.ID);
Labels = vertcat(S.Labels);
Y2 = Labels == 0;
Y1 = ~Y2;
Y = Y1 - Y2;
mean_Y = mean(Y);
X = X - repmat(mean(X), size(X,1), 1);
Ym = Y - mean(Y);

[b_acc_opls,n_LV,w_opls,yhat_opls] = DAMACY_top(X,ID,Y,5,0);
% w(ind) = w;

%% Plot
if plotme == 1
    VariableNames = horzcat(VariableNames,VariableNames);
    Analysis_Result_1D(S_A,w_opls,edges,binsize,b_acc_opls,yhat_opls,Y,VariableNames,'OPLS',dataname);
end

%% Top Model SVM
% SVM_model = svmtrain(Labels,X,'-c 700 -t 2 -g 2 -r 1 -q');
% [yhatt_svm, b_acc_svm, yhat_svm] = svmpredict(Labels, X, SVM_model, '');
% [b_acc] = DAMACY_SVM_top(X,ID,Labels);
% w_svm = SVM_model.SVs' * SVM_model.sv_coef;
% Analysis_Result_1D(w_svm,edges,binsize,b_acc_svm(1)/100,yhat_svm,VariableNames, 'SVM',dataname);
end



