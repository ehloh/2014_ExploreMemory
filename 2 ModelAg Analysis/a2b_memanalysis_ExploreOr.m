% Is memory better when subjects choose to explore? 
clear all; clc; close all hidden;

% Choose
log.taskversion='Both';
log.SecondHalfOnly=0; 
%
log.v1_1_ok={'p01_AF' ; 'p01_AS' ; 'p02_TY' ; 'p03_JV' ; 'p06_MK' ; 'p08_AK' ; 'p10_SG' ; 'p11_OG' ; 'p30_AR' ; 'p31_RD' ; 'p33_AA' ; 'p34_BE'};
log.v1_2_ok={'p14_CP';'p15_SZ';'p23_SG';'p24_VP';'p26_PB';'p27_PV';'p29_MY';'p35_EO';'p37_ND';'p38_SL';'p39_JH';'p40_CC'; 'p45_SO';};
log.specificsubjects=log.v1_1_ok;
% log.specificsubjects=[log.v1_1_ok; log.v1_2_ok];

for o1=1:1 % General setup
    
    % Folders
    w=pwd; if strcmp(w(2), ':')==1; where.where='D:\Dropbox\SCRIPPS\8 Explore Mem'; else; where.where='/Users/EleanorL/Dropbox/SCRIPPS/8 Explore Mem';  end
    where.data=[where.where filesep '2 All data' filesep log.taskversion]; % Append version
    where.analysisscript=[where.where filesep '3a Analysis scripts'];
    where.analysisinputs=[where.where filesep '4 Inputs'];
    if isdir(where.analysisinputs)==0; mkdir(where.analysisinputs); end
    addpath(where.analysisscript)
    
    % Fetch subjects
    [w.s w.s1 log.koshertable]=xlsread([where.analysisscript filesep 'i_subjectsok.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.koshertable, log.specificsubjects, log.koshertable, 'Include');
    
    % Mark  condition
    for s=1:log.n_subjs
        if sum(strcmp(log.subjects{s},log.v1_1_ok))==1; log.subj_condition(s,1)=1;
        elseif sum(strcmp(log.subjects{s},log.v1_2_ok))==1; log.subj_condition(s,1)=2;
        else error(['What condition is this subject in ??? ' log.subjects{s}]); 
        end
    end
    
    disp(' ##################################################')
    disp(['No. of subjects:  ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;  disp(log.subjects); end
    input('Hit enter to continue     '); disp(' ')
    disp(' ##################################################')
end

%% Fetch data + calculate mem FA rates 

% subjdata: (1) name (2) raw data & data in cells (3) data split into cells (4) memory scores per cell
    % Col 1=Name, Col 2=Raw data & data in cells, Col 3=Data split into cells, Col 4=memory scores epr cell
subjdata=cell(log.n_subjs,2); %d_subjstats=cell(log.n_subjs,2);
for s=1:log.n_subjs
    subjdata{s,1}=log.subjects{s};
    
    % Load encoding data
    ws.enc=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_3ExploreEnc.mat']);
    subjdata{s,2}.col.e=ws.enc.cF.col;
    subjdata{s,2}.encdata=ws.enc.cF.rdata;
    subjdata{s,2}.design=ws.enc.cF.par.design;
    
    % Load memory data
    ws.mem=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_4Memtest.mat']);
    subjdata{s,2}.memdata=ws.mem.mem.data;
    subjdata{s,2}.col=ws.mem.mem.col;
    
    % Collate data (memdata)
    if s==1;  % columns
        col.e=ws.enc.cF.col; 
        col=ws.mem.mem.col;  % 'col' refers to memtest trialstats
        %
        col.NTokens=2;
        col.pLoss=24;
        col.Entropy=25;
        col.EntropyNTok=26;
        col.EV=27;
        col.OutcomeMagnitude=col.OutcomeMag;
        col.OutcomeMean=28;
        col.OutcomeVariance=29;
        col.NTokenPairs=[];
        col.ChoiceExplore=30; % Binary
        
        % Insert new column/criteria here **
        
        
    end
    subjdata{s,2}.memdata(:, col.NTokens)  =subjdata{s,2}.memdata(:, col.NTokens)*2;
    ws.new=subjdata{s,2}.memdata( subjdata{s,2}.memdata(:, col.Item_OldNew)==2,:) ;
    ws.new(:, [col.NTokens col.pLoss col.Entropy col.EntropyNTok col.EV col.OutcomeMean col.OutcomeVariance col.ChoiceExplore])=nan;
    ws.old=subjdata{s,2}.memdata( subjdata{s,2}.memdata(:, col.Item_OldNew)==1,:);
    ws.old=ws.old(ws.old(:,col.Enc_TypeOK)==1,:); % correctly encoded trials only
    ws.old(:, col.ChoiceExplore)=0; ws.old(ws.old(:, col.Choice)==3, col.ChoiceExplore)=1;
    
    % Write parameters into memdata
    [ ws.old ] = fpar_conflict( ws.old , col);
    subjdata{s,2}.memdata= sortrows([ws.old; ws.new],  col.Mem_Trialnum);
    subjdata{s,2}.olds=ws.old;
    subjdata{s,2}.foils=ws.new;
    try subjdata{s,2}.behstats=ws.enc.cF.behstats; catch, disp(['Stats missing for ' log.subjects{s}]); end
    
    % Memory FA rates 
    subjdata{s,3}.FA.hit=sum(subjdata{s,2}.foils(:,col.OldNew)==1)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.surehit=sum(subjdata{s,2}.foils(:,col.OldNew)==1 & subjdata{s,2}.foils(:,col.SureGuess)==1)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.rem=sum(subjdata{s,2}.foils(:,col.RemKnow)==1)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.know=sum(subjdata{s,2}.foils(:,col.RemKnow)==2)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.surerem=sum(subjdata{s,2}.foils(:,col.RemKnow)==1 & subjdata{s,2}.foils(:,col.SureGuess)==1)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.sureknow=sum(subjdata{s,2}.foils(:,col.RemKnow)==2 & subjdata{s,2}.foils(:,col.SureGuess)==1)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.Pos1=sum(subjdata{s,2}.foils(:,col.Pos)==1)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.Pos2=sum(subjdata{s,2}.foils(:,col.Pos)==2)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.Pos3=sum(subjdata{s,2}.foils(:,col.Pos)==3)/size(subjdata{s,2}.foils,1);
    subjdata{s,3}.FA.Pos4=sum(subjdata{s,2}.foils(:,col.Pos)==4)/size(subjdata{s,2}.foils,1);
    
    %
    ws=[];
end


% Overall scores (sanity check)
d_meanmem=nan(log.n_subjs, 1); 
for s=1:log.n_subjs
    k=1; 
    ws.mold=  subjdata{s,2}.memdata(subjdata{s,2}.memdata(:, col.Item_OldNew)==1,:);  
    
    % Hit
    d_meanmem(s,k)= mean( ws.mold(:, col.OldNew) ==1 ) ;     k=k+1; 
    
    % Sure Hit 
    d_meanmem(s,k)= mean( ws.mold(:, col.OldNew)==1 & ws.mold(:, col.SureGuess)==1 ) ;       k=k+1; 
    
    % Rem  
    d_meanmem(s,k)= mean( ws.mold(:, col.RemKnow)==1 ) ;     k=k+1; 
    
    % Know 
    d_meanmem(s,k)= mean( ws.mold(:, col.RemKnow)==2  ) ;     k=k+1; 
    
    % Sure Rem  
    d_meanmem(s,k)= mean( ws.mold(:, col.RemKnow)==1 & ws.mold(:, col.SureGuess)==1 ) ;     k=k+1; 
    
    % Sure Know 
    d_meanmem(s,k)= mean( ws.mold(:, col.RemKnow)==2 & ws.mold(:, col.SureGuess)==1 ) ;     k=k+1; 
    
    
    % Assoc mem
    ws.correctrec= ws.mold(ws.mold(:, col.CorrectRecognition)==1,:); 
    d_meanmem(s,k)= mean(ws.correctrec(:, col.CorrectPosition));       k=k+1;   
    
    
end
% 

mean(d_meanmem)
std(d_meanmem)

error



%% Calculate memory scores

% Which cells to include for ExploreOr analysis?
for o=1:1
    
    % Instruc cells: Subject, In-cluster cells (env threat, ntokens)
%     instruc.exploreor_groupcells=[2 3; 2 4;3 3]; % 0.4< p(Explore) <0.6  (mean)
    instruc.exploreor_groupcells=[1 4; 2 3; 2 4; 3 2; 3 3; 3 4]; % 0.25 < p(Explore) <0.6  (mean)  
    
%     instruc.exploreor_groupcells=[ 1 1; 1 2; 1 3; 1 4; 2 1; 2 2; 2 3; 2 4; 3 1; 3 2; 3 3; 3 4; 4 1; 4 2; 4 3; 4 4;]; % ALL CELLS
    
    
    %
    instruc.exploreor_cells=[log.subjects repmat({instruc.exploreor_groupcells}, log.n_subjs,1)];
    disp('Cells included in ExploreOr sample (EnvThreat, NTok):');  disp(instruc.exploreor_groupcells)

end

% Sort data in cells + calculate memory scores (according to instruc.cells)
%       instruc.exploreor_cells: col 3 = n Explore trials, col 4 = n Or trials 
for o1=1:1 % Results tables 't_': 1st col=Explore, 2nd col=Or
    log.subjectsexcluded={}; % Checks 
    w.res=nan*zeros(log.n_subjs,2);  % Results tables 
    t_dprime=w.res;
    t_hit =w.res;
    t_surehit =w.res;
    t_rem =w.res;
    t_know =w.res;
    t_surerem =w.res;
    t_sureknow =w.res;
    t_assoc =w.res;
    t_phitsure =w.res;
    t_phitrem =w.res;
    t_phitsurerem =w.res;
    t_premsure =w.res;
    t_pknowsure =w.res;
end
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ': ' log.subjects{s}])
    ws.d=subjdata{s,2}.olds;
    
    % Which cells for this subject?
    ws.whichcells=instruc.exploreor_cells{find(strcmp(instruc.exploreor_cells(:,1), log.subjects{s})),2};
        
    % Are there instructions for this subject? Or, on valid trials?
    if isempty(ws.whichcells)==0  
       
        % Automated compile data, including specified cells 
        ws.sample=[];
        for c=1:size(ws.whichcells,1)
            wc.d= ws.d(ws.d(:,col.EnvThreat)*4 ==ws.whichcells(c,1) &  ws.d(:,col.NTokens)/2==ws.whichcells(c,2) , :) ;
            ws.sample=[ws.sample; wc.d];
            wc=[];
        end
        ws.Explore=ws.sample(ws.sample(:,col.Choice)==3,:);
        ws.Or=ws.sample(ws.sample(:,col.Choice)~=3,:);
        subjdata{s,3}.ExploreOrSample=ws.sample;
        subjdata{s,3}.ExploreTrials=ws.Explore;
        subjdata{s,3}.OrTrials=ws.Or;
        instruc.exploreor_cells{find(strcmp(instruc.exploreor_cells(:,1), log.subjects{s})),3}=size(ws.Explore,1);
        instruc.exploreor_cells{find(strcmp(instruc.exploreor_cells(:,1), log.subjects{s})),4}=size(ws.Or,1);
        
        % Log results ---------------------------
        % (1) dprime
        ws.mem=f_calcmem(1, ws.Explore(:,col.Roc), subjdata{s,3}.FA.hit);
        t_dprime(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(1, ws.Or(:,col.Roc), subjdata{s,3}.FA.hit);
        t_dprime(s,2)=ws.mem.rate;
        % (2) rem
        ws.mem=f_calcmem(2, ws.Explore(:,col.Roc), subjdata{s,3}.FA.rem);
        t_rem(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(2, ws.Or(:,col.Roc), subjdata{s,3}.FA.rem);
        t_rem(s,2)=ws.mem.rate;
        % (3) know
        ws.mem=f_calcmem(3, ws.Explore(:,col.Roc), subjdata{s,3}.FA.know);
        t_know(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(3, ws.Or(:,col.Roc), subjdata{s,3}.FA.know);
        t_know(s,2)=ws.mem.rate;
        % (4) surerem
        ws.mem=f_calcmem(4, ws.Explore(:,col.Roc), subjdata{s,3}.FA.surerem);
        t_surerem(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(4, ws.Or(:,col.Roc), subjdata{s,3}.FA.surerem);
        t_surerem(s,2)=ws.mem.rate;
        % (5) sureknow
        ws.mem=f_calcmem(5, ws.Explore(:,col.Roc), subjdata{s,3}.FA.sureknow);
        t_sureknow(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(5, ws.Or(:,col.Roc), subjdata{s,3}.FA.sureknow);
        t_sureknow(s,2)=ws.mem.rate;
        % (6) hit
        ws.mem=f_calcmem(6, ws.Explore(:,col.Roc), subjdata{s,3}.FA.hit);
        t_hit(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(6, ws.Or(:,col.Roc), subjdata{s,3}.FA.hit);
        t_hit(s,2)=ws.mem.rate;
        % (7) surehit
        ws.mem=f_calcmem(7, ws.Explore(:,col.Roc), subjdata{s,3}.FA.surehit);
        t_surehit(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(7, ws.Or(:,col.Roc), subjdata{s,3}.FA.surehit);
        t_surehit(s,2)=ws.mem.rate;
        % (8) Associative memory
        t_assoc(s,1)=sum( ws.Explore(:,col.CorrectPosition) ) / sum( ws.Explore(:, col.OldNew)==1 );
        t_assoc(s,2)=sum( ws.Or(:,col.CorrectPosition) ) / sum( ws.Or(:, col.OldNew)==1 );
        % (9) percent Hit that is sure
        ws.mem=f_calcmem(8, ws.Explore(:,col.Roc),0);
        t_phitsure(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(8, ws.Or(:,col.Roc),0);
        t_phitsure(s,2)=ws.mem.rate;
        % (10) percent Hit that is Rem
        ws.mem=f_calcmem(9, ws.Explore(:,col.Roc),0);
        t_phitrem(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(9, ws.Or(:,col.Roc),0);
        t_phitrem(s,2)=ws.mem.rate;
        % (11) percent Hit that is Surerem
        ws.mem=f_calcmem(10, ws.Explore(:,col.Roc),0);
        t_phitsurerem(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(10, ws.Or(:,col.Roc),0);
        t_phitsurerem(s,2)=ws.mem.rate;
        % (12) percent Rem that is Sure
        ws.mem=f_calcmem(11, ws.Explore(:,col.Roc),0);
        t_premsure(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(11, ws.Or(:,col.Roc),0);
        t_premsure(s,2)=ws.mem.rate;
        % (13) percent Know that is Sure
        ws.mem=f_calcmem(12, ws.Explore(:,col.Roc),0);
        t_pknowsure(s,1)=ws.mem.rate;
        ws.mem=f_calcmem(12, ws.Or(:,col.Roc),0);
        t_pknowsure(s,2)=ws.mem.rate;
        
        %
        ws=[];
    elseif isempty(ws.whichcells)==0 
        disp('No trials for this subject. If unexpected, check excel input')
        log.subjectsexcluded{size(log.subjectsexcluded,1)+1,1}=log.subjects{s};
    else error('Could not find instructions for this subject, specifying which cells are Explore and which are Or')
    end
    
end
disp('Check no. of Explore/Or memory trials  [  Col 3-4= Explore & Or, no. trials  ]'); disp(instruc.exploreor_cells)
% disp(['Mean no. trials (' num2str( size(instruc.exploreor_groupcells,1) )  ' cells) : Explore = ' num2str(mean(cell2mat(instruc.exploreor_cells(:,3))),3) ', Or = '  num2str(mean(cell2mat(instruc.exploreor_cells(:,4))),3)])
disp(['Mean no. trials (' num2str( size(instruc.exploreor_groupcells,1) )  ' cells) : Explore = ' num2str(mean(cell2mat(instruc.exploreor_cells(:,3))),3)  '(sd=' num2str(std(cell2mat(instruc.exploreor_cells(:,3))),3)  '), Or = '  num2str(mean(cell2mat(instruc.exploreor_cells(:,4))),3) ' (sd=' num2str(std(cell2mat(instruc.exploreor_cells(:,4))),3) ')'])

mean(cell2mat(instruc.exploreor_cells(:,3:4)))

std(cell2mat(instruc.exploreor_cells(:,3:4)))

%% Does choice correlate w stuff in the limited task space?

checkvars={'EnvThreat';'NTokens';'pLoss';'Entropy';'EntropyNTok';'EV'; 'OutcomeMagnitude';};


% Frank correlations
d_varcorexp=[[{'Subject'} checkvars']; [log.subjects cell(log.n_subjs, length(checkvars))]];
for c=1:length(checkvars)
    eval([ 'wv.colnum=col.' checkvars{c} ';']);
    for s=1:log.n_subjs
        [r p]=corr(subjdata{s,3}.ExploreOrSample(:, [wv.colnum col.ChoiceExplore]));
        [r p]=corr(subjdata{s,3}.ExploreOrSample(:, [wv.colnum col.ChoiceExplore]), 'type', 'Kendall');
        
        if p<0.1
        d_varcorexp{s+1, c+1}=r;
        else d_varcorexp{s+1, c+1}='-';
        end
    end

end

% GLM
dothis=0;
if dothis
    stopped here. the best way to check this is really to use the behavioural modelling. so, dont do this 
    d_varcorexp=[[{'Subject'} checkvars']; [log.subjects cell(log.n_subjs, length(checkvars))]];
    for c=1:length(checkvars)
        eval([ 'wv.colnum=col.' checkvars{c} ';']);
        for s=1:log.n_subjs
            [r p]=corr(subjdata{s,3}.ExploreOrSample(:, [wv.colnum col.ChoiceExplore]));
            [r p]=corr(subjdata{s,3}.ExploreOrSample(:, [wv.colnum col.ChoiceExplore]), 'type', 'Kendall');
            
            if p<0.1
                d_varcorexp{s+1, c+1}=r;
            else d_varcorexp{s+1, c+1}='-';
            end
        end
        
    end
end

%% Statistical test (within condition)

instruc.memtypes={'dprime'; 'hit';'surehit'; 'assoc';'rem';'know';'surerem';'sureknow';};
instruc.memtypes_names={'Recognition memory (d'')'; 'Hit';'Sure hit'; 'Position recall accuracy';'Remember';'Know';'Sure remember';'Sure know';};
%         'phitrem';'phitsure';'phitsurerem';'premsure';'pknowsure'};
r_exploreor=[instruc.memtypes num2cell(zeros(length(instruc.memtypes),3))];  % h p tstat
for m=1:length(instruc.memtypes)
    eval(['wm.mem=t_' instruc.memtypes{m} ';']);
    [h p wm.ci wm.stats]=ttest(wm.mem(:,1)-wm.mem(:,2));
    r_exploreor(m, 2:3)=[{h} {p}]; 
    if p<0.1;  r_exploreor{m,2}=0.5; end
    
    r_exploreor(m, 2:end)=[{h} {p} {wm.stats.tstat}]; % Print full details
end
% openvar r_exploreor

% % plot some only
% instruc.memtypes={'dprime'; 'assoc';'rem';'know';}; instruc.memtypes_names={'Recognition memory (d'')'; 'Position recall accuracy';'Remember';'Know'};
instruc.memtypes={'dprime';  'rem';'know';}; instruc.memtypes_names={'Recognition memory (d'')'; 'Remember rates';'Know rates'};

% 

% Plot
close all hidden
f.plotcols=2;  f.figwidth= 1200; f.figheight=500; f.fontsize=45; f.fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms, Cambria
f.subplot_VerHorz=[0.1 0.1]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.1 0.1];
figure('Name', 'ExploreOr memory', 'NumberTitle', 'off', 'Position', [100 50 f.figwidth f.figheight], 'Color', 'w');  k=1;
for m=1:length(instruc.memtypes)
    eval(['wm.mem=t_' instruc.memtypes{m} ';']);
    
    subplot(ceil(length(instruc.memtypes)/f.plotcols),f.plotcols, k)
    barwitherr(std(wm.mem)/sqrt(log.n_subjs),  mean(wm.mem),'m'  )  % Mean +/- Std Error
    title(instruc.memtypes_names{m}, 'FontSize', f.fontsize, 'FontName', f.fontname)
    set(gca, 'xtick', 1:2,'xticklabel', {'Explore';'Or'}, 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
    set(gca, 'xtick', 1:2,'xticklabel', {'Explore';'Accept/Reject'}, 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
    xlim([0.5 2.5])
    k=k+1;
end


%% Compare conditions

% Compile data
instruc.printmem={'dprime' 'd'; 'hit' 'h';'surehit' 'sh'; 'assoc' 'a';
    'rem' 'r';'know' 'k';'surerem' 'sr';'sureknow' 'sk';};
%     'phitrem' 'phr';'phitsure' 'phs';'phitsurerem' 'phsr';'premsure' 'prs';'pknowsure' 'pks'};
d_memscores=[{'Subject' 'Condition'}; [log.subjects num2cell(log.subj_condition)]];
for m=1:size(instruc.printmem,1)
    eval(['wm.mem=t_' instruc.printmem{m,1} ';'])
    d_memscores=[d_memscores [{[instruc.printmem{m,2} '_Explore'] [instruc.printmem{m,2} '_Or']}; num2cell(wm.mem)]];
end
for r=2:size(d_memscores,1) % Remove nans
    for c=2:size(d_memscores,2)
        if isnan(d_memscores{r,c})==1
            d_memscores{r,c}=[];
        end
    end
end
