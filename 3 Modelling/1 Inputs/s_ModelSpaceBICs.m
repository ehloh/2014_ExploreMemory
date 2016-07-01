% Look at BICs over the entire model space
clear all; close all hidden; clc

cd('2b Hierar Fits');
load('res_hierarfitmodels_cF (16-Apr-2015) v1 bpjm10_L3066.mat')
% load('res_hierarfitmodels_cF (16-Apr-2015) v2 bpj16_L3489.mat')

% Compare w Bayes factor?
m1=1;
m2=2;
B=exp((r_res{m1,3}-r_res{m2,3})*-0.5);
disp(['B=' num2str(B)  '  (m1=' r_res{m1,1} ', m2='  r_res{m2,1} ')     [3-10=moderate evidence, >10=strong]'])

%%
r_res=sortrows(r_res,-3);
bics=cell2mat(r_res(:,3));
n_params=cellfun(@(x)size(x,2)-3, r_res(:,2));
modnames=r_res(:,1);

% Sorting by # params
all=[num2cell(bics) num2cell(n_params) modnames]; all=sortrows(all, 2);
trans=[1 find(cell2mat(all(:,2))==2,1,'first') find(cell2mat(all(:,2))==3,1,'first') find(cell2mat(all(:,2))==4,1,'first') find(cell2mat(all(:,2))==5,1,'first') find(cell2mat(all(:,2))==6,1,'first') find(cell2mat(all(:,2))==7,1,'first')];

% Plot
% figure('color','w'); k=1; f.FontSize=20; f.FontName='PT Sans';
% subplot(4,1,k); bar(bics-min(bics)); title('Difference in BIC from winning model'); k=k+1; xlabel('Model'), ylabel('Change in BIC')
% subplot(4,1,k); bar(bics./n_params); title('BIC/No. params'); k=k+1; xlabel('Model'), ylabel('BIC/No. Params')
% subplot(4,1,k); 





figure('color', 'w'), f.FontSize=25; f.FontName='PT Sans';
bar(cell2mat(all(:,1)) ); title('BIC (ranked by no. params)', 'FontSize', f.FontSize, 'FontName', f.FontName)
k=k+1; xlabel('Model (grouped by no. of parameters)', 'FontSize', f.FontSize, 'FontName', f.FontName)
ylabel('BIC', 'FontSize', f.FontSize, 'FontName', f.FontName)
set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
for i=2:7
    hx = graph2d.constantline(trans(i)-0.5, 'Color',[0 0 0]); changedependvar(hx,'x');
end
hold on, bar( find(cell2mat(all(:,1))== min(cell2mat(all(:,1)))), min(bics), 'r'); % highlist best 
ylim([3000 4800])
set(gca, 'xtick', [])








