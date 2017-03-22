function Analysis_Result_1D(S,w,edges,binsize,b_acc,yhat,Y,VariableNames,regmeth,dataname,paired)
%% Created by Thomas Vijverberg on 09-06-2015 at Radboud University Nijmegen
% Last edited by Thomas Vijverberg on 09-06-2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input data structure (STRUCT) 'S' with dimension 1 x (total_measurements)
    %Must contain atleast contain:
    %Data input (S.Data)
    %IDs input (S.ID)
    %Labels S.Labels

%Input weights vector OPLS (vect) 1x(var*(binsize-2))

%Input matrix edges from histogram (MAT)(var x binsize-2)

%Input binsize for histograms (INT) 20-1500 recommended

%Input accuracy OPLS (double)

%Input yhat OPLS vector (vect) (length(Y) x 1)

%Input Labels vector (vect) (length(Y) x 1)

%Input VariableNames matrix (MAT) with dimension 1 x (number_variables)

%Input string regression method (STR)

%Input string data name (STR)

%Input if data is paired (BOOL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
N_var = (length(w)/(binsize-2));
N_ID = length(S);
ID = vertcat(S.ID);
top_cutoff = 0.015;
bot_cutoff = 0.01;
if paired == 1
    paired = 'LOOCV';
else
    paired = 'LTOCV';
end

load('cmap_bluered.mat');

%% Vectorize weights
for i = 1:(length(w)/(binsize-2))
    vector_end = [1+(i-1)*(binsize-2) i*(binsize-2)];
    w_var(i,:) = w(vector_end(1):vector_end(2));
end

max_w = max(w_var');
min_w = min(w_var');

for i = 1:N_var
    if max_w(i) < top_cutoff
        top_cutoff = max_w(i)-0.005;
    end
    if min_w(i) > -1*bot_cutoff
        bot_cutoff = -1*min_w(i)+0.005;
    end     
end


%% Edgeline
w_edge_high = w_var;
w_edge_low = w_var;
w_edge_high2 = w_var;
w_edge_low2 = w_var;

for i = 1:N_var
    w_edge_high2(i,:) = w_var(i,:) > max_w(i)-top_cutoff;
    w_edge_low2(i,:) = w_var(i,:) < min_w(i)+bot_cutoff;
    w_edge_high(i,:) = w_edge_high2(i,:).*w_var(i,:);
    w_edge_low(i,:) = w_edge_low2(i,:).*w_var(i,:);
end


%% Area Under Curve (AUC)
for i = 1:N_var
    auc(i) = trapz(abs(w_var(i,:)));
end

tot_auc = sum(auc);
auc = auc ./tot_auc;

for i = 1:N_var
    auc_str{i} = [' ' num2str(round(auc(i)*100)) '%'];
end
auc_pie_str = strcat(VariableNames,auc_str);

%% Plotting
% --------------------Variable weight plots----------------------
figure
hold on
for i = 1:N_var
    subtightplot(N_var,2,i*2-1,[0.05 0.1]);
    plot(edges(i,2:binsize-1),w_var(i,:),'g');
    line(edges(i,2:binsize-1),repmat(max_w(i)-top_cutoff,binsize-2,1));
    line(edges(i,2:binsize-1),repmat(min_w(i)+bot_cutoff,binsize-2,1));
    texts = ['Variable ' VariableNames{i}];
    title(texts,'Fontsize' ,8);
    ylim([min(min_w) max(max_w)]);
    if i == round(length(edges(:,1))/2)
        ylabel(['Marker weights of ' dataname])
    end
end
xlabel('Autoscaled Marker Intensity');

% --------------------AUC Pie plot----------------------
subtightplot(1,4,4);
pie(auc,auc_pie_str);
title('Marker influence (AUC)', 'Fontsize',16);
hold on;

idx1 = Y == -1;
idx2 = Y == 1;
yhat = yhat - mean([mean(yhat(idx2)) mean(yhat(idx1))]);
subplot(1, 4 ,3, 'Ylim', [-2 2]);
hold on;
plot(0,yhat(idx1),'r*')
text(repmat(0.25,sum(idx1),1), yhat(idx1), num2str(ID(1:sum(idx1))), 'VerticalAlignment','top','HorizontalAlignment','left', 'Fontsize', 5, 'color','r')
plot(0,yhat(idx2), 'b*')
text(repmat(0.25,sum(idx2),1), yhat(idx2), num2str(ID(sum(idx1)+1:length(Y))), 'VerticalAlignment','top','HorizontalAlignment','left', 'Fontsize', 5, 'color','b')
line([-1 1], [0 0])
title('yhat of controls (red) vs diseased (blue)', 'fontsize', 22)

axes('Units','Normal');
h = title(['DAMACY of ' dataname ' is predicted with ' num2str(b_acc*100) '% accuracy (' regmeth ', ' paired ')'], 'fontsize', 22 ,'Position',[0.5 1.05 0.5]);
set(gca,'visible','off')
set(h,'visible','on')

% --------------------Original Variable Plot ----------------------
figure
hold on

data = vertcat(S.Data);
w_combin = w_edge_high2+(-1*w_edge_low2);
for i = 1:N_var
    his(i) = histogram(data(:,i),binsize);
    his_max(i) = max(his(i).Values);
end
his_m = max(his_max);

for i = 1:N_var
    subtightplot(N_var,1,i,[0.05 0.1]);
    hold on;
    his = histogram(data(:,i),binsize);
    h=imagesc([edges(i,1) edges(i,end)], [0 his_m+100],w_combin(i,:));
    alpha(h,0.3);
    ylim([0 his_m+100]);
    colormap(cmap);
    texts = ['Variable ' VariableNames{i}];
    title(texts,'Fontsize' ,8);
    if i == round(length(edges(:,1))/2)
        ylabel(['Number of Cells of ' dataname])
    end
end
xlabel('Autoscaled Marker Intensity');
axes('Units','Normal');
h = title(['DAMACY histograms & phenotype regions of ' dataname '. Red is control; blue is case, ' num2str(b_acc*100) '% accuracy (' regmeth ', ' paired ')'], 'fontsize', 22 ,'Position',[0.5 1.05 0.5]);
set(gca,'visible','off')
set(h,'visible','on')

% --------------------Contour plot----------------------
%for l1 = 1:length(S)
% figure;
% count = 1;
% for i = 1:N_var
%     for j = 1:N_var
%         if j > i && mod(count-1,7)
%             subtightplot(N_var,N_var,count,[0.04 0.03]);
%             Histogram(1,:,:,:,:,:,:) = NDhist(S(1).Data(:,[i j]), edges([i j],:), 3);
%             colormap(flipud(gray));
%             nom = imagesc(edges(i,:),edges(j,:),squeeze(Histogram(1,:,:,:,:,:,:,:))');
%             alpha(nom,0.5);
%             hold on;
%             freezeColors();
%             
%             I1 = repmat(w_edge_high(i,:),binsize-2,1);
%             J1 =repmat(w_edge_high(j,:),binsize-2,1);
%             IJ1 = I1 + J1';
%             h1 = imagesc(edges(i,2:end-1),edges(j,2:end-1),IJ1');
%             colormap(cmap(33:64,:));
%             alpha(h1,0.05);
%             freezeColors();
%             
%             I2 = repmat(w_edge_low(i,:),binsize-2,1);
%             J2 =repmat(w_edge_low(j,:),binsize-2,1);
%             IJ2 = I2 + J2';
%             h2 = imagesc(edges(i,2:end-1),edges(j,2:end-1),IJ2');
%             colormap(flipud(cmap(1:32,:)));
%             alpha(h2,0.05);
%             xlabel(VariableNames{i});
%             ylabel(VariableNames{j});
%             %                high = w_edge_high(i,:)'*w_edge_high(j,:);
%             %                low = w_edge_low(i,:)'*w_edge_low(j,:);
%             %                contour(edges(i,2:end-1),edges(j,2:end-1),high,'b');
%             %                contour(edges(i,2:end-1),edges(j,2:end-1),low,'r');
%             xlabel(VariableNames{i});
%             ylabel(VariableNames{j});
%         end
%         count = count +1;
%     end
% end
end

