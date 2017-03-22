%GRIDSEARCH_HISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

dataname = input('What is the name of the dataset?\n','s');
grid = [80];

for j = 1:8
    St = Get_Tube(S,j);
    for i = 1:length(grid)
        [b_acc_opls(j,i)] = Complete_Histograms( St,grid(i),VariableNames,dataname,0);
        disp(['Assesment of gridsize ' num2str(grid(i)) ' resulted in: ' num2str(b_acc_opls(j,i)*100) '% accuracy in Tube ' num2str(j)]);
    end
    disp('');
end

[b_acc,ind] = max(b_acc_opls);
[b_acc,tube] = max(b_acc);

%[b_acc]=Complete_Histograms(S,grid(ind),VariableNames,dataname,1)
