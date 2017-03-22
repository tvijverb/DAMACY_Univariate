function [ N ] = Create_histogram( S,var )
%CREATE_HISTOGRAM multidimensional histogram
%   S is datainput [X,var] matrix
%   var is amount of variables to use in histogram (integer)
load('50-grid.mat');
N = zeros(length(S),length(edges)-1);
temp = zeros(length(S),length(edges)-1);

%Filling first variable into histogram matrix
for i = 1:length(S);
    N(i,:) = histcounts(S(i).Data(:,1),edges);
end

%Adding following variables in 3D space
if var >= 2
    for j = 2:var
        for i = 1:length(S)
            temp(i,:) = histcounts(S(i).Data(:,j),edges);  
        end
     N = cat(3,N,temp);
end
end

