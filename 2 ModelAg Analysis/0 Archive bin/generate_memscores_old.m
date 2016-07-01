 d_a_surerem=w.grid; d_a_know=w.grid; d_a_sureknow=w.grid;
d_a_phitsure=w.grid; d_a_phitrem=w.grid; d_a_phitsurerem=w.grid; d_a_premsure=w.grid; d_a_pknowsure=w.grid;
d_r_dprime=w.grid; d_r_hit=w.grid; d_r_surehit=w.grid; d_r_assoc=w.grid;
d_r_rem=w.grid; d_r_surerem=w.grid; d_r_know=w.grid; d_r_sureknow=w.grid;
d_r_phitsure=w.grid; d_r_phitrem=w.grid; d_r_phitsurerem=w.grid; d_r_premsure=w.grid; d_r_pknowsure=w.grid;
d_e_dprime=w.grid; d_e_hit=w.grid; d_e_surehit=w.grid; d_e_assoc=w.grid;
d_e_rem=w.grid; d_e_surerem=w.grid; d_e_know=w.grid; d_e_sureknow=w.grid;


% Sort data in cells + calculate memory scores (according to instruc.cellstions)
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ': ' log.subjects{s}])
    ws.d=subjdata{s,2}.olds;
    subjdata{s,3}.celldata=cell(size(instruc.cells,1),2); % data sorted into cells (col 1=name, col 2= trialstats)
    
    % Automated sort into cells
    for i=1:size(instruc.cells,1)
        disp(['Cell: ' instruc.cells{i,1}])
        
        % Compile sample
         wc.d=ws.d;
        for c=1:size(instruc.cells{i,2},1) % apply each criteria in series
            disp(['      Criteria ' num2str(c) ':  ' instruc.cells{i,2}{c,1} ' == '  num2str(instruc.cells{i,2}{c,2})])
            eval(['wc.d=wc.d(wc.d(:, col.m.' instruc.cells{i,2}{c,1} ')==' num2str(instruc.cells{i,2}{c,2}) ',:);']);
        end
        disp(['     Total no items in cell: ' num2str(length(wc.d(:,col.m.Roc)))])
        subjdata{s,3}.celldata{i,1}=instruc.cells{i,1};
        subjdata{s,3}.celldata{i,2}=wc.d;
       
        % Results tables ---------------------------        
        % (1) dprime
        wc.mem=f_calcmem(1, wc.d(:,col.m.Roc), subjdata{s,3}.FA.hit);
        t_dprime{s+1,i+1}=wc.mem.rate;
           if isnan(wc.mem.rate);
            t_dprime{s+1,i+1}=w.empty;
        else
        t_dprime{s+1,i+1}=wc.mem.rate;        
           end
        % (2) rem
        wc.mem=f_calcmem(2, wc.d(:,col.m.Roc), subjdata{s,3}.FA.rem);
        t_rem{s+1,i+1}=wc.mem.rate;
        if isnan(wc.mem.rate);
            t_rem{s+1,i+1}=w.empty;
        else
        t_rem{s+1,i+1}=wc.mem.rate;        
           end
        % (3) know
        wc.mem=f_calcmem(3, wc.d(:,col.m.Roc), subjdata{s,3}.FA.know);
        t_know{s+1,i+1}=wc.mem.rate;
           if isnan(wc.mem.rate);
            t_know{s+1,i+1}=w.empty;
        else
        t_know{s+1,i+1}=wc.mem.rate;        
           end
        % (4) surerem
        wc.mem=f_calcmem(4, wc.d(:,col.m.Roc), subjdata{s,3}.FA.surerem);
        t_surerem{s+1,i+1}=wc.mem.rate;
          if isnan(wc.mem.rate);
            t_surerem{s+1,i+1}=w.empty;
        else
        t_surerem{s+1,i+1}=wc.mem.rate;        
          end
        % (5) sureknow
        wc.mem=f_calcmem(5, wc.d(:,col.m.Roc), subjdata{s,3}.FA.sureknow);
        t_sureknow{s+1,i+1}=wc.mem.rate;
           if isnan(wc.mem.rate);
            t_sureknow{s+1,i+1}=w.empty;
        else
        t_sureknow{s+1,i+1}=wc.mem.rate;        
           end
        % (6) hit
        wc.mem=f_calcmem(6, wc.d(:,col.m.Roc), subjdata{s,3}.FA.hit);
        t_hit{s+1,i+1}=wc.mem.rate;
            if isnan(wc.mem.rate);
            t_hit{s+1,i+1}=w.empty;
        else
        t_hit{s+1,i+1}=wc.mem.rate;        
           end
        % (7) surehit
        wc.mem=f_calcmem(7, wc.d(:,col.m.Roc), subjdata{s,3}.FA.surehit);
        t_surehit{s+1,i+1}=wc.mem.rate;
            if isnan(wc.mem.rate);
            t_surehit{s+1,i+1}=w.empty;
        else
        t_surehit{s+1,i+1}=wc.mem.rate; 
            end
        % (8) Associative memory
        t_assoc{s+1,i+1}=sum(wc.d(:,col.m.CorrectPosition))/size(wc.d,1);
        % (9) percent Hit that is sure
        wc.mem=f_calcmem(8, wc.d(:,col.m.Roc),0);
        t_phitsure{s+1,i+1}=wc.mem.rate;
        % (10) percent Hit that is Rem
        wc.mem=f_calcmem(9, wc.d(:,col.m.Roc),0);
        t_phitrem{s+1,i+1}=wc.mem.rate;
        % (11) percent Hit that is Surerem
        wc.mem=f_calcmem(10, wc.d(:,col.m.Roc),0);
        t_phitsurerem{s+1,i+1}=wc.mem.rate;
        % (12) percent Rem that is Sure
        wc.mem=f_calcmem(11, wc.d(:,col.m.Roc),0);
        t_premsure{s+1,i+1}=wc.mem.rate;
        % (13) percent Know that is Sure
        wc.mem=f_calcmem(12, wc.d(:,col.m.Roc),0);
        t_pknowsure{s+1,i+1}=wc.mem.rate;
        
        % Data grids (if cell involves EnvThreat x NTokens classification)--------------------------------
        if sum(strcmp(instruc.cells{i,2}(:,1), 'EnvThreat_Level'))~=0 & sum(strcmp(instruc.cells{i,2}(:,1), 'NTokens_Level'))~=0
            e=length(design.NTokenPairs_Npairs)+1- instruc.cells{i,2}{find(strcmp(instruc.cells{i,2}(:,1),'EnvThreat_Level')),2}; % grid allocation
            n=instruc.cells{i,2}{find(strcmp(instruc.cells{i,2}(:,1),'NTokens_Level')),2};
            
            % Write memory scores to grid
            wc.dprime=t_dprime{s+1,i+1};
            wc.hit=t_hit{s+1,i+1};
            wc.surehit=t_surehit{s+1,i+1};
            wc.rem=t_rem{s+1,i+1};
            wc.surerem=t_surerem{s+1,i+1};
            wc.know=t_know{s+1,i+1};
            wc.sureknow=t_sureknow{s+1,i+1};
            wc.assoc=t_assoc{s+1,i+1};
            wc.phitsure=t_phitsure{s+1,i+1};
            wc.phitrem=t_phitrem{s+1,i+1};
            wc.phitsurerem=t_phitsurerem{s+1,i+1};
            wc.premsure=t_premsure{s+1,i+1};
            wc.pknowsure=t_pknowsure{s+1,i+1};
            
            % Choice is included in cell criteria
            
            
            if sum(strcmp(instruc.cells{i,2}(:,1), 'Choice'))~=0
                if instruc.cells{i,2}{find(strcmp(instruc.cells{i,2}(:,1), 'Choice')),2}==1 % Accept
                    d_a_dprime{e,n}{s}=wc.dprime;
                    d_a_hit{e,n}{s}=wc.hit;
                    d_a_surehit{e,n}{s}=wc.surehit;
                    d_a_rem{e,n}{s}=wc.rem;
                    d_a_surerem{e,n}{s}=wc.surerem;
                    d_a_know{e,n}{s}=wc.know;
                    d_a_sureknow{e,n}{s}=wc.sureknow;
                    d_a_assoc{e,n}{s}=wc.assoc;
                    d_a_phitsure{e,n}{s}=wc.phitsure;
                    d_a_phitrem{e,n}{s}=wc.phitrem;
                    d_a_phitsurerem{e,n}{s}=wc.phitsurerem;
                    d_a_premsure{e,n}{s}=wc.premsure;
                    d_a_pknowsure{e,n}{s}=wc.pknowsure;
                elseif instruc.cells{i,2}{find(strcmp(instruc.cells{i,2}(:,1), 'Choice')),2}==2 % Reject
                    d_r_dprime{e,n}{s}=wc.dprime;
                    d_r_hit{e,n}{s}=wc.hit;
                    d_r_surehit{e,n}{s}=wc.surehit;
                    d_r_rem{e,n}{s}=wc.rem;
                    d_r_surerem{e,n}{s}=wc.surerem;
                    d_r_know{e,n}{s}=wc.know;
                    d_r_sureknow{e,n}{s}=wc.sureknow;
                    d_r_assoc{e,n}{s}=wc.assoc;
                    d_r_phitsure{e,n}{s}=wc.phitsure;
                    d_r_phitrem{e,n}{s}=wc.phitrem;
                    d_r_phitsurerem{e,n}{s}=wc.phitsurerem;
                    d_r_premsure{e,n}{s}=wc.premsure;
                    d_r_pknowsure{e,n}{s}=wc.pknowsure;
                elseif instruc.cells{i,2}{find(strcmp(instruc.cells{i,2}(:,1), 'Choice')),2}==3 % Explore
                    d_e_dprime{e,n}{s}=wc.dprime;
                    d_e_hit{e,n}{s}=wc.hit;
                    d_e_surehit{e,n}{s}=wc.surehit;
                    d_e_rem{e,n}{s}=wc.rem;
                    d_e_surerem{e,n}{s}=wc.surerem;
                    d_e_know{e,n}{s}=wc.know;
                    d_e_sureknow{e,n}{s}=wc.sureknow;
                    d_e_assoc{e,n}{s}=wc.assoc;
                    d_e_phitsure{e,n}{s}=wc.phitsure;
                    d_e_phitrem{e,n}{s}=wc.phitrem;
                    d_e_phitsurerem{e,n}{s}=wc.phitsurerem;
                    d_e_premsure{e,n}{s}=wc.premsure;
                    d_e_pknowsure{e,n}{s}=wc.pknowsure;
                end
            else % Independent of choice
                d_dprime{e,n}{s}=wc.dprime;
                d_hit{e,n}{s}=wc.hit;
                d_surehit{e,n}{s}=wc.surehit;
                d_rem{e,n}{s}=wc.rem;
                d_surerem{e,n}{s}=wc.surerem;
                d_know{e,n}{s}=wc.know;
                d_sureknow{e,n}{s}=wc.sureknow;
                d_assoc{e,n}{s}=wc.assoc;
                d_phitsure{e,n}{s}=wc.phitsure;
                d_phitrem{e,n}{s}=wc.phitrem;
                d_phitsurerem{e,n}{s}=wc.phitsurerem;
                d_premsure{e,n}{s}=wc.premsure;
                d_pknowsure{e,n}{s}=wc.pknowsure;
            end
        end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        wc=[];
    end
    
    ws=[];
end

%% Output memory scores to txt


if request.printmemscores
    instruc.printmem={'dprime';'rem';'know';'surerem';'sureknow';'hit';'surehit'; 'assoc';'phitrem';'phitsure';'phitsurerem';'premsure';'pknowsure'};
    for m=1:length(instruc.printmem)
        eval(['wm.d=t_' instruc.printmem{m} ';'])
        
        
        for r=2:size(wm.d,1)
            for c=2:size(wm.d,2)
                if isnan(wm.d{r,c})==1
                    wm.d{r,c}='';
                end
            end
        end
        
        
        [printok]=print2txt(where.analysisinputs, ['(' date ') ' instruc.printmem{m}], wm.d);
        disp([instruc.printmem{m} ' ------- '])
        disp(printok)
    end
    
end

%% Grids: Calculate mean & plot

% Which memory measure
instruc.gridmem={
    'dprime';'hit';'surehit';'rem';'surerem';'know';'sureknow';'assoc' ; 'phitsure';'phitrem';'phitsurerem';'premsure';'pknowsure';
    'a_dprime';'a_hit';'a_surehit';'a_rem';'a_surerem';'a_know';'a_sureknow';'a_assoc'; %'a_phitsure';'a_phitrem';'a_phitsurerem';'a_premsure';'a_pknowsure';
    'r_dprime';'r_hit';'r_surehit';'r_rem';'r_surerem';'r_know';'r_sureknow';'r_assoc'; %'r_phitsure';'r_phitrem';'r_phitsurerem';'r_premsure';'r_pknowsure';
    'e_dprime';'e_hit';'e_surehit';'e_rem';'e_surerem';'e_know';'e_sureknow';'e_assoc'; %'e_phitsure';'e_phitrem';'e_phitsurerem';'e_premsure';'e_pknowsure';
    };

% Calculate means
for i=1:length(instruc.gridmem)
    eval(['wg.d=d_' instruc.gridmem{i} ';'])
    wg.r=cell(size(wg.d));
    
    
    
    for e=1:size(wg.d,1)
        for n=1:size(wg.d,2)
            if sum(isnan(cell2mat(wg.d{e,n} ))==0)>=log.grid_minimumsubjs % minimum no. subjects?
                wg.r{size(wg.d,1)+1-e,n}=mean( cell2mat(wg.d{e,n}(isnan(cell2mat(wg.d{e,n} ))==0)) );
                wg.sd{size(wg.d,1)+1-e,n}=std(cell2mat(wg.d{e,n}(isnan(cell2mat(wg.d{e,n} ))==0)));
            else
                wg.r{size(wg.d,1)+1-e,n}=nan;
                wg.sd{size(wg.d,1)+1-e,n}=nan;
            end
        end
    end
    
    eval(['m_' instruc.gridmem{i} '=wg.r;'])
    eval(['sd_' instruc.gridmem{i} '=wg.sd;'])
    wg=[];
end


% Plot (specify)
which='dprime';
%
eval(['m_mem= m_' which ';'])
eval(['sd_mem= sd_' which ';'])
%
close all hidden; figure; set(gcf,'Color',[1 1 1])
design.EnvThreat_labels=cell(length(design.EnvThreat_pBomb),1); % labels
for i=1:length(design.EnvThreat_pBomb) % Hard-coded!!
    
%     
%     [n d]=rat(design.EnvThreat_pBomb)
%     ;
    
    [n d]=rat(design.EnvThreat_pBomb(i));
    design.EnvThreat_labels{length(design.EnvThreat_pBomb)+1-i}=[ num2str(n) '/' num2str(d)];
end
subplot(1,2,1); imagesc(cell2mat(m_mem));     % MEAN
colorbar; title(['mean ' which]); axis square
ylabel('EnvThreat');  set(gca,'YTick',1:length(design.EnvThreat_labels)); set(gca, 'YTickLabel', design.EnvThreat_labels)
xlabel('N Tokens'); set(gca,'XTick',1:length(design.NTokenPairs_Npairs)); set(gca, 'XTickLabel',design.NTokenPairs_Npairs)
subplot(1,2,2); imagesc(cell2mat(sd_mem));        % SD
colorbar; title(['sd ' which]); axis square
ylabel('EnvThreat'); set(gca,'YTick',1:length(design.EnvThreat_labels)); set(gca, 'YTickLabel', design.EnvThreat_labels)
xlabel('N Tokens'); set(gca,'XTick',1:length(design.NTokenPairs_Npairs)); set(gca, 'XTickLabel',design.NTokenPairs_Npairs)
disp(['Means for displayed mem score ' which]); disp(m_mem); disp(' ')


% LEFT TO DO:
% - Abigail: write script to combine cells (use eval functions)
% - Manually plot all grids, and present

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               % Generate scores for all cells in design
clear all; clc; close all hidden;

% where.where='E:\Analysis';
where.where='/Volumes/NO NAME/Analysis';

log.taskversion='v1_1';
log.specificsubjects={'p01_AF' ; 'p01_AS' ; 'p02_TY' ; 'p03_JV' ; 'p06_MK' ; 'p08_AK' ; 'p10_SG' ; 'p11_OG' ; 'p30_AR' ; 'p31_RD' ; 'p33_AA' ; 'p34_BE'};
%log.specificsubjects={'p14_CP' ; 'p15_SZ' ; 'p23_SG' ; 'p24_VP' ; 'p26_PB' ; 'p27_PV' ; 'p29_MY' ; 'p35_EO' ; 'p37_ND' ; 'p38_SL' ; 'p39_JH' ; 'p40_CC'};
log.grid_minimumsubjs=2;

request.printmemscores=0;

for o1=1:1 % General setup
    
    % Folders
    where.data=[where.where filesep '2 All data' filesep log.taskversion]; % Append version
    where.analysisscript=[where.where filesep '3 Analysis scripts'];
    where.analysisinputs=[where.where filesep '4 Inputs'];
    if isdir(where.analysisinputs)==0; mkdir(where.analysisinputs); end
    addpath(where.analysisscript)
    
    % Fetch subjects
 
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep '3 Analysis scripts' filesep 'i_subjectsok.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.koshertable, log.specificsubjects, log.koshertable, 'Include');
    
%   cd(where.data); log.subjects=dir('p*'); log.subjects={log.subjects(:).name}';
%    for i=1:length(log.specificsubjects)
%        if sum(sum(strcmp(log.subjects, log.specificsubjects{i}))==1)
%            log.subjects{find(strcmp(log.subjects, log.specificsubjects{i})),2}=1;
%        else
%            error(['Unidentified specific subject requested: ' log.specificsubjects{i}]);
%        end
%     end
%    if isempty(log.specificsubjects)==1; log.subjects(:,2)=num2cell(ones(size(log.subjects,1),1)); end
%    log.subjects=log.subjects(cell2mat(log.subjects(:,2))==1,1);
%    log.n_subjs=length(log.subjects);
    
    % Instructions
    %   Col 1:                  Name of cell
    %   Col 2 onwards:    {Factor name, Value}
    log.taskversion(log.taskversion=='-')='_';
    [w.a instruc.cells w.b]=xlsread([where.analysisscript filesep 'memanalysis_instruc.xlsx'], log.taskversion);
    for i=1:size(instruc.cells,1)
        eval(['instruc.cells{i,2}=' instruc.cells{i,2} ';'])
    end

    disp(' ##################################################')
    disp(['No. of subjects:  ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0
        disp(log.subjects)
    else disp(log.specificsubjects)
    end
    input('Hit enter to continue     '); disp(' ')
    disp(' ##################################################')
end

w.empty=nan;

%% Fetch data (encoding and mem)

% Fetch data from all subjects, and task parameters from 1st subject (col, design) 
    % Col 1=Name, Col 2=Raw data & data in cells, Col 3=Data split into cells, Col 4=memory scores epr cell
subjdata=cell(log.n_subjs,2);
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
    subjdata{s,2}.col.m=ws.mem.mem.col;
    
    % Collate data (memdata)
    if s==1;  % columns
        col.e=ws.enc.cF.col; 
        col.m=ws.mem.mem.col; 
        design=ws.enc.cF.par.design;
        %
        col.m.EnvThreat_Level=21;
        col.m.NTokens_Level=22;
        col.m.NTokens=23;
        col.m.pLoss=24;
        col.m.Entropy=25;
        col.m.EVExplore=26;
        col.m.EV=27;
        col.m.OutcomeMean=28;
        col.m.OutcomeVariance=29;
        
        % Insert new column/criteria here **
        
    end
    ws.new=subjdata{s,2}.memdata( subjdata{s,2}.memdata(:, col.m.Item_OldNew)==2,:) ;
    ws.new(:, [col.m.NTokens col.m.pLoss col.m.Entropy col.m.EVExplore col.m.EV col.m.OutcomeMean col.m.OutcomeVariance])=nan;
    ws.old=subjdata{s,2}.memdata( subjdata{s,2}.memdata(:, col.m.Item_OldNew)==1,:) ;
    ws.old(:,col.m.EnvThreat_Level)=ws.old(:,col.m.EnvThreat);
    for i=1:length( subjdata{s,2}.design.NTokenPairs_Npairs)
        ws.old(ws.old(:,col.m.NTokenPairs)==i, col.m.NTokens_Level)=i;
    end
    ws.old=ws.old(ws.old(:,col.m.Enc_TypeOK)==1,:); % correctly encoded trials only
    
    % Write parameters into memdata
    [ ws.old ws.col] = fpar_conflict( ws.old , col.m, ws.enc.cF.par.design);
    %
    subjdata{s,2}.memdata= sortrows([ws.old; ws.new],  col.m.Mem_Trialnum);
    subjdata{s,2}.olds=ws.old;
    subjdata{s,2}.foils=ws.new;
    
    %
    ws=[];
end


%% Mark additional variables (on subjdata{s,2}.olds)

design.var.pLoss=unique(subjdata{1,2}.olds(:,col.m.pLoss));
design.var.Entropy=unique(subjdata{1,2}.olds(:,col.m.Entropy));
design.var.EVExplore=unique(subjdata{1,2}.olds(:,col.m.EVExplore));
design.var.EV=unique(subjdata{1,2}.olds(:,col.m.EV));

col.m.pLoss_Third=30;
col.m.pLoss_Third1=31;
col.m.pLoss_Third1_Accept=32;
col.m.pLoss_Third1_Reject=33;
col.m.pLoss_Third1_Explore=34;
col.m.pLoss_Third2=35;
col.m.pLoss_Third2_Accept=36;
col.m.pLoss_Third2_Reject=37;
col.m.pLoss_Third2_Explore=38;
col.m.pLoss_Third3=39;
col.m.pLoss_Third3_Accept=40;
col.m.pLoss_Third3_Reject=41;
col.m.pLoss_Third3_Explore=42;
col.m.pLoss_Med=43;
col.m.pLoss_Med1=44;
col.m.pLoss_Med1_Accept=45;
col.m.pLoss_Med1_Reject=46;
col.m.pLoss_Med1_Explore=47;
col.m.pLoss_Med2=48;
col.m.pLoss_Med2_Accept=49;
col.m.pLoss_Med2_Reject=50;
col.m.pLoss_Med2_Explore=51;
col.m.InCluster=52;
col.m.Choice_Explore=53;
col.m.InCluster_Explore=54;
col.m.InCluster_NonExplore=55;
col.m.EnvThreat_Level_Med=56;
col.m.EnvThreat_Level_Med1=57;
col.m.EnvThreat_Level_Med1_Accept=58;
col.m.EnvThreat_Level_Med1_Reject=59;
col.m.EnvThreat_Level_Med1_Explore=60;
col.m.EnvThreat_Level_Med2=61;
col.m.EnvThreat_Level_Med2_Accept=62;
col.m.EnvThreat_Level_Med2_Reject=63;
col.m.EnvThreat_Level_Med2_Explore=64;
col.m.NTokens_Level_Med=65;
col.m.NTokens_Level_Med1=66;
col.m.NTokens_Level_Med1_Accept=67;
col.m.NTokens_Level_Med1_Reject=68;
col.m.NTokens_Level_Med1_Explore=69;
col.m.NTokens_Level_Med2=70;
col.m.NTokens_Level_Med2_Accept=71;
col.m.NTokens_Level_Med2_Reject=72;
col.m.NTokens_Level_Med2_Explore=73;
col.m.NTokens_Med=74;
col.m.NTokens_Med1=75;
col.m.NTokens_Med1_Accept=76;
col.m.NTokens_Med1_Reject=77;
col.m.NTokens_Med1_Explore=78;
col.m.NTokens_Med2=79;
col.m.NTokens_Med2_Accept=80;
col.m.NTokens_Med2_Reject=81;
col.m.NTokens_Med2_Explore=82;
col.m.Entropy_Med=83;
col.m.Entropy_Med1=84;
col.m.Entropy_Med1_Accept=85;
col.m.Entropy_Med1_Reject=86;
col.m.Entropy_Med1_Explore=87;
col.m.Entropy_Med2=88;
col.m.Entropy_Med2_Accept=89;
col.m.Entropy_Med2_Reject=90;
col.m.Entropy_Med2_Explore=91;
col.m.Entropy_Quart=92;
col.m.Entropy_Quart1=93;
col.m.Entropy_Quart1_Accept=94;
col.m.Entropy_Quart1_Reject=95;
col.m.Entropy_Quart1_Explore=96;
col.m.Entropy_Quart2=97;
col.m.Entropy_Quart2_Accept=98;
col.m.Entropy_Quart2_Reject=99;
col.m.Entropy_Quart2_Explore=100;
col.m.Entropy_Quart3=101;
col.m.Entropy_Quart3_Accept=102;
col.m.Entropy_Quart3_Reject=103;
col.m.Entropy_Quart3_Explore=104;
col.m.Entropy_Quart4=105;
col.m.Entropy_Quart4_Accept=106;
col.m.Entropy_Quart4_Reject=107;
col.m.Entropy_Quart4_Explore=108;
col.m.EVExplore_Med=109;
col.m.EVExplore_Med1=110;
col.m.EVExplore_Med1_Accept=111;
col.m.EVExplore_Med1_Reject=112;
col.m.EVExplore_Med1_Explore=113;
col.m.EVExplore_Med2=114;
col.m.EVExplore_Med2_Accept=115;
col.m.EVExplore_Med2_Reject=116;
col.m.EVExplore_Med2_Explore=117;
col.m.EVExplore_Third=118;
col.m.EVExplore_Third1=119;
col.m.EVExplore_Third1_Accept=120;
col.m.EVExplore_Third1_Reject=121;
col.m.EVExplore_Third1_Explore=122;
col.m.EVExplore_Third2=123;
col.m.EVExplore_Third2_Accept=124;
col.m.EVExplore_Third2_Reject=125;
col.m.EVExplore_Third2_Explore=126;
col.m.EVExplore_Third3=127;
col.m.EVExplore_Third3_Accept=128;
col.m.EVExplore_Third3_Reject=129;
col.m.EVExplore_Third3_Explore=130;
col.m.EV_Med=131;
col.m.EV_Med1=132;
col.m.EV_Med1_Accept=133;
col.m.EV_Med1_Reject=134;
col.m.EV_Med1_Explore=135;
col.m.EV_Med2=136;
col.m.EV_Med2_Accept=137;
col.m.EV_Med2_Reject=138;
col.m.EV_Med2_Explore=139;
col.m.EV_Third=140;
col.m.EV_Third1=141;
col.m.EV_Third1_Accept=142;
col.m.EV_Third2_Reject=143;
col.m.EV_Third3_Explore=144;
col.m.EV_Third2=145;
col.m.EV_Third2_Accept=146;
col.m.EV_Third2_Reject=147;
col.m.EV_Third2_Explore=148;
col.m.EV_Third3=149;
col.m.EV_Third3_Accept=150;
col.m.EV_Third3_Reject=151;
col.m.EV_Third3_Explore=152;
col.m.Dummy=153;




for s=1:log.n_subjs
    ws.d=subjdata{s,2}.olds;
    
    for t=1:size(ws.d,1)
          
        % Learning variables
        for o1=1:1 % pLoss
            % pLoss Third
            if ws.d(t,col.m.pLoss)<=design.var.pLoss(3)
                ws.d(t,col.m.pLoss_Third)=1;
            elseif ws.d(t,col.m.pLoss)<=design.var.pLoss(6) && ws.d(t,col.m.pLoss)>design.var.pLoss(3)
                ws.d(t,col.m.pLoss_Third)=2;
            elseif ws.d(t,col.m.pLoss)>design.var.pLoss(6)
                ws.d(t,col.m.pLoss_Third)=3;
            else
                error('Unidentified pLoss value!')
            end
            
            % pLoss Median
            if ws.d(t,col.m.pLoss)<=design.var.pLoss(5)
                ws.d(t,col.m.pLoss_Med)=1;
            elseif ws.d(t,col.m.pLoss)>design.var.pLoss(5)
                ws.d(t,col.m.pLoss_Med)=2;
            else
                error('Unidentified pLoss value!')
            end
        end
        
        for o2=1:1 % EnvThreat_Level
             % EnvThreat_Level Median
            if ws.d(t,col.m.EnvThreat_Level)<=2
                ws.d(t,col.m.EnvThreat_Level_Med)=1;
            elseif ws.d(t,col.m.EnvThreat_Level)>2
                ws.d(t,col.m.EnvThreat_Level_Med)=2;
            else error('Unidentified EnvThreat_Level value!')
            end
            
        end
        
        for o3=1:1 % NTokens_Level    
            % NTokens_Level Median
            if ws.d(t,col.m.NTokens_Level)<=2
                ws.d(t,col.m.NTokens_Level_Med)=1;
            elseif ws.d(t,col.m.NTokens_Level)>2
                ws.d(t,col.m.NTokens_Level_Med)=2;
            else error('Unidentified NTokens_Level value!')
            end
            
        end    
        
        for o4=1:1 % NTokens
            % NToken Median
            if ws.d(t,col.m.NTokens)<=4
                ws.d(t,col.m.NTokens_Med)=1;
            elseif ws.d(t,col.m.NTokens)>4
                ws.d(t,col.m.NTokens_Med)=2;
            else error('Unidentified NTokens value!')
            end

        end
        
        for o5=1:1 % Entropy
            % Entropy Median
            if ws.d(t,col.m.Entropy)<=design.var.Entropy(4)
                ws.d(t,col.m.Entropy_Med)=1;
            elseif ws.d(t,col.m.Entropy)>design.var.Entropy(4)
                ws.d(t,col.m.Entropy_Med)=2;
            else
                error('Unidentified Entropy value!')
            end
            
            % Entropy Quarters
            if ws.d(t,col.m.Entropy)<=design.var.Entropy(2)
                ws.d(t,col.m.Entropy_Quart)=1;
            elseif ws.d(t,col.m.Entropy)<=design.var.Entropy(4) && ws.d(t,col.m.Entropy)>design.var.Entropy(2)
                ws.d(t,col.m.Entropy_Quart)=2;
            elseif ws.d(t,col.m.Entropy)<=design.var.Entropy(6) && ws.d(t,col.m.Entropy)>design.var.Entropy(4)
                ws.d(t,col.m.Entropy_Quart)=3;
            elseif ws.d(t,col.m.Entropy)>design.var.Entropy(6)
                ws.d(t,col.m.Entropy_Quart)=4;
            else
                error('Unidentified Entropy value!')
            end
            
        end
        
        for o6=1:1 % EVExplore
            % EVExplore Median
            if ws.d(t,col.m.EVExplore)<=design.var.EVExplore(8)
                ws.d(t,col.m.EVExplore_Med)=1;
            elseif ws.d(t,col.m.EVExplore)>design.var.EVExplore(8)
                ws.d(t,col.m.EVExplore_Med)=2;
            else
                error('Unidentified EVExplore value!')
            end
            
            % EVExplore Third
             if ws.d(t,col.m.EVExplore)<=design.var.EVExplore(5)
                ws.d(t,col.m.EVExplo