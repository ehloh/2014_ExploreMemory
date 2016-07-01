% Plot individual & group choice
clear all; close all hidden; clc

% where.where='/Volumes/SANDISK/8 Explore Mem';
where.where='D:\Dropbox\SANDISK\8 Explore Mem';

%%
log.v1_1_ok={'p01_AF' ; 'p01_AS' ; 'p02_TY' ; 'p03_JV' ; 'p06_MK' ; 'p08_AK' ; 'p10_SG' ; 'p11_OG' ; 'p30_AR' ; 'p31_RD' ; 'p33_AA' ; 'p34_BE'};
log.v1_2_ok={'p14_CP';'p15_SZ';'p23_SG';'p24_VP';'p26_PB';'p27_PV';'p29_MY';'p35_EO';'p37_ND';'p38_SL';'p39_JH';'p40_CC'; 'p45_SO';};


log.TaskVersion='Both';
log.SecondHalfOnly=0;

% 
% log.specificsubjects={};
log.specificsubjects=log.v1_2_ok;
% log.specificsubjects=[log.v1_1_ok ; log.v1_2_ok];

for o1=1:1 % General setup  
    
    where.alldata=[where.where filesep '2 All Data'];
    where.execscript=([where.where filesep '1 Explore Mem Task']); addpath(where.execscript)
    addpath([where.where filesep '3a Analysis scripts'])
    
    % Fetch subjects
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep '3a Analysis scripts' filesep 'i_subjectsok.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.koshertable, log.specificsubjects, log.koshertable, 'Include');
       
end

%% Load all subject data (and append details)

disp('Loading subject data ----------------------------------')
subjdata=cell(log.n_subjs,3); % Col 1= Name, Col 2= Trialstats, Col 3=Choice grids, Row n+1 = Group
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ' (' log.subjects{s} ') -- '])
    ws=load([where.alldata filesep log.TaskVersion filesep log.subjects{s} filesep log.subjects{s} '_file_3ExploreEnc.mat']);
    if s==1; col=ws.cF.col; design=ws.cF.par.design; end
    
    % Append details to Encoding file (behstats)
    ws.cF.behstats.per_a=sum(ws.cF.rdata(ws.cF.rdata(:,col.Trialvalid)==1,col.Resp1)==1)/sum(ws.cF.rdata(:,col.Trialvalid));
    ws.cF.behstats.per_r=sum(ws.cF.rdata(ws.cF.rdata(:,col.Trialvalid)==1,col.Resp1)==2)/sum(ws.cF.rdata(:,col.Trialvalid));
    ws.cF.behstats.per_e=sum(ws.cF.rdata(ws.cF.rdata(:,col.Trialvalid)==1,col.Resp1)==3)/sum(ws.cF.rdata(:,col.Trialvalid));
    ws.cF.behstats.totalwin= sum(ws.cF.rdata(:,col.OutcomeMag));
    cF=[];cF=ws.cF;
    save([where.alldata filesep log.TaskVersion filesep log.subjects{s} filesep log.subjects{s} '_file_3ExploreEnc.mat'], 'cF');
    
    % Second half only?
    if log.SecondHalfOnly
        ws.cF.rdata=ws.cF.rdata(round(size(ws.cF.rdata,1)/2)+1:end,:);
        [ws.cF.behaviour] = f_plotindividualbehaviour(ws.cF.par.design,ws.cF.rdata, col); close; 
    end
    
    
    
    % Record data
    subjdata{s,1}=log.subjects{s};
    subjdata{s,2}=ws.cF.rdata;
    subjdata{s,3}.a=ws.cF.behaviour.accept;
    subjdata{s,3}.r=ws.cF.behaviour.reject;
    subjdata{s,3}.e=ws.cF.behaviour.explore;
    
    % Mark group
    if sum(strcmp( log.subjects{s},log.v1_1_ok)); subjdata{s,4}=1;
    elseif sum(strcmp( log.subjects{s},log.v1_2_ok)); subjdata{s,4}=2;
    end
        
        %
    ws=[];
end



% Table for output to SPSS
t_accept=[[{'Subject'}; log.subjects] [{'Condition'};   subjdata(:,4)]];  t_reject=t_accept;t_explore=t_accept;kk=3;
for s=1:log.n_subjs
    k=kk;


for e=1:4
    ee=5-e;
    for n=1:4
        if s==1; t_accept{1,k}=['a_e' num2str(e) 'n' num2str(n)]; t_reject{1,k}=['r_e' num2str(e) 'n' num2str(n)]; t_explore{1,k}=['e_e' num2str(e) 'n' num2str(n)];end
        
        
        t_accept{s+1,k}=subjdata{s,3}.a(ee,n);
        t_reject{s+1,k}=subjdata{s,3}.r(ee,n);
        t_explore{s+1,k}=subjdata{s,3}.e(ee,n);
        
        
        
        k=k+1;
    end
end

end
t_choice=[t_accept t_reject(:, 3:end) t_explore(:, 3:end)];

% Create mean choice stats (average choice stats across subjects)
subjdata{log.n_subjs+1,3}.a=zeros(design.EnvThreat_Levels, design.NTokenPairs_Levels);
subjdata{log.n_subjs+1,3}.r=zeros(design.EnvThreat_Levels, design.NTokenPairs_Levels);
subjdata{log.n_subjs+1,3}.e=zeros(design.EnvThreat_Levels, design.NTokenPairs_Levels);
for s=1:log.n_subjs
    subjdata{log.n_subjs+1,3}.a=subjdata{log.n_subjs+1,3}.a+subjdata{s,3}.a;
    subjdata{log.n_subjs+1,3}.r=subjdata{log.n_subjs+1,3}.r+subjdata{s,3}.r;
    subjdata{log.n_subjs+1,3}.e=subjdata{log.n_subjs+1,3}.e+subjdata{s,3}.e;
end
subjdata{log.n_subjs+1,3}.a=subjdata{log.n_subjs+1,3}.a/log.n_subjs;
subjdata{log.n_subjs+1,3}.r=subjdata{log.n_subjs+1,3}.r/log.n_subjs;
subjdata{log.n_subjs+1,3}.e=subjdata{log.n_subjs+1,3}.e/log.n_subjs;
subjdata{log.n_subjs+1,1}='Group mean';

%% Plot p(Choice) 

disp('Subjects included ################'); disp(char(log.subjects))

% Generate labels
for e=1:design.EnvThreat_Levels
    f.ylabels{design.EnvThreat_Levels+1-e}=strtrim(rats(design.EnvThreat_pBomb(e)));
end
for n=1:design.NTokenPairs_Levels
    f.xlabels{n}=num2str(2*design.NTokenPairs_Npairs(n));
end


fontsize=15; fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms


for o1=1:1  % Plot group 
    close all hidden
    figure('Name', ['Group Choice plots (' log.TaskVersion ', n= ' num2str(log.n_subjs) ')'],'Position',[200,00, 1200, 600], 'NumberTitle','off'); set(gcf,'Color', 'w');
    
    subplot(1,3,1) % Accept
    imagesc(subjdata{log.n_subjs+1,3}.a, [0 1])
    axis('square'); title('% Accept','FontSize',15, 'FontName', fontname); colorbar
    ylabel('EnvThreat'); set(gca,'YTick',1:design.EnvThreat_Levels)
    set(gca, 'YTickLabel', f.ylabels)
    xlabel('NTokens');
    set(gca,'XTick',1:design.NTokenPairs_Levels)
    set(gca, 'XTickLabel',f.xlabels)
    ylabel('Environmental Threat','FontSize',15, 'FontName', fontname)
    xlabel('No. Tokens','FontSize',15, 'FontName', fontname)
    set(gca,'FontSize',15, 'FontName', fontname)
    
    subplot(1,3,2) % Reject
    imagesc(subjdata{log.n_subjs+1,3}.r, [0 1])
    axis('square'); title('% Reject','FontSize',15, 'FontName', fontname); colorbar
    ylabel('EnvThreat');
    ylabel('EnvThreat'); set(gca,'YTick',1:design.EnvThreat_Levels)
    set(gca, 'YTickLabel', f.ylabels)
    xlabel('NTokens');
    set(gca,'XTick',1:design.NTokenPairs_Levels)
    set(gca, 'XTickLabel',f.xlabels)
    ylabel('Environmental Threat','FontSize',15, 'FontName', fontname)
    xlabel('No. Tokens','FontSize',15, 'FontName', fontname)
    set(gca,'FontSize',15, 'FontName', fontname)
    
    subplot(1,3,3) % Explore
    imagesc(subjdata{log.n_subjs+1,3}.e, [0 1])
    axis('square'); title('% Explore','FontSize',15, 'FontName', fontname); colorbar
    ylabel('EnvThreat');
    ylabel('EnvThreat'); set(gca,'YTick',1:design.EnvThreat_Levels)
    set(gca, 'YTickLabel', f.ylabels)
    xlabel('NTokens');
    set(gca,'XTick',1:design.NTokenPairs_Levels)
    set(gca, 'XTickLabel',f.xlabels)
    ylabel('Environmental Threat','FontSize',15, 'FontName', fontname)
    xlabel('No. Tokens','FontSize',15, 'FontName', fontname)
    set(gca,'FontSize',15, 'FontName', fontname)
end


% Plot individuals 
f.figheight=log.n_subjs*80; 
f.subplotcols=4;   f.subplot_VerHorz=[0.0095 0.07]; f.fig_BotTop=[0.005 0.025]; f.submw=[0.05 0.1]; 
f.printdetails_xy=[0.8 1.2];
f.figindex=  figure('Name', ['Individual Choice plots (' log.TaskVersion ', n= ' num2str(log.n_subjs) ')'],'Position',[200,50, 900, f.figheight], 'NumberTitle','off'); set(gcf,'Color', 'w');
for s=1:log.n_subjs % Plot subjects
   
    % Subject name
    subtightplot(log.n_subjs, f.subplotcols,  (s-1)*f.subplotcols+1, f.subplot_VerHorz, f.fig_BotTop, f.submw); 
    text(1.1,0.4, log.subjects{s}); axis off
    if s==1; text(f.printdetails_xy(1), f.printdetails_xy(2), ['Task version: ' log.TaskVersion]); axis off; end
    
    % Accept
    subtightplot(log.n_subjs, f.subplotcols,  (s-1)*f.subplotcols+2, f.subplot_VerHorz, f.fig_BotTop, f.submw); 
    imagesc(subjdata{s,3}.a, [0 1])
    axis('square'); axis off; if s==1; title('% Accept'); end
    
    % Reject
    subtightplot(log.n_subjs, f.subplotcols,  (s-1)*f.subplotcols+3, f.subplot_VerHorz, f.fig_BotTop, f.submw); 
    imagesc(subjdata{s,3}.r, [0 1])
    axis('square'); axis off; if s==1; title('% Reject'); end
    
    % Explore
    subtightplot(log.n_subjs, f.subplotcols,  (s-1)*f.subplotcols+4, f.subplot_VerHorz, f.fig_BotTop, f.submw); 
    imagesc(subjdata{s,3}.e, [0 1])
    axis('square'); axis off; if s==1; title('% Explore'); end; colorbar
end

%%

d_meanchoice=[cellfun(@(x) mean(x(:, col.Resp1)==1), subjdata(1: log.n_subjs,2))  cellfun(@(x) mean(x(:, col.Resp1)==2), subjdata(1: log.n_subjs,2)) cellfun(@(x) mean(x(:, col.Resp1)==3), subjdata(1: log.n_subjs,2))];
d_mc1=d_meanchoice(find(cellfun(@(x)sum(strcmp(log.v1_1_ok, x)), log.subjects)), :);
d_mc2=d_meanchoice(find(1-cellfun(@(x)sum(strcmp(log.v1_1_ok, x)), log.subjects)), :);


for c=1:3
[h p]=ttest2(d_mc1(:,c), d_mc2(:,c)); 
disp([ num2str(c) '  ' num2str(p)]);
end


stds=[std(d_mc1(:,1))/ sqrt(size(d_mc1,1))  std(d_mc2(:,1))/ sqrt(size(d_mc2,1))
    std(d_mc1(:,2))/ sqrt(size(d_mc1,1))  std(d_mc2(:,2))/ sqrt(size(d_mc2,1));
    std(d_mc1(:,3))/ sqrt(size(d_mc1,1))  std(d_mc2(:,3))/ sqrt(size(d_mc2,1))]; 
means=[mean(d_mc1(:,1)) mean(d_mc2(:,1))
    mean(d_mc1(:,2)) mean(d_mc2(:,2))
    mean(d_mc1(:,3)) mean(d_mc2(:,3))];


figure('color','w')
barwitherr(stds, means)
set(gca, 'FontSize', fontsize, 'FontName', fontname); 
legend({'Immediate test';'Delayed test'})
set(gca, 'xticklabel', {'% Accept','% Reject','% Explore'}, 'FontSize', fontsize, 'FontName', fontname); 
ylabel( 'Percentage choice', 'FontSize', fontsize, 'FontName', fontname); 
title('Overall proportion of choices', 'FontSize', fontsize, 'FontName', fontname); 




%%

close all hidden


ans=num2str(subjdata{log.n_subjs+1,3}.e,3);



