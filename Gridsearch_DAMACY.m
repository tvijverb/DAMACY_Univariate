%GRIDSEARCH_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here
dataname = input('What is the name of the dataset?\n','s');
grid = [20 30 40 50 80 150 250 400 600 800 1000];

for i = 1:length(grid)
    [b_acc_opls(i)] = Complete_Histograms( S,grid(i),VariableNames,dataname,0);
    disp(['Assesment of gridsize ' num2str(grid(i)) ' resulted in: ' num2str(b_acc_opls(i)*100) '% accuracy']);
end

[b_acc,ind] = max(b_acc_opls);

[b_acc]=Complete_Histograms(S,grid(ind),VariableNames,dataname,1)
