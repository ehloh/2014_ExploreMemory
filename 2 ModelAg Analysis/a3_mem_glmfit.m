% Generate scores for all cells in design
% clear all; clc; close all hidden;
clear all; clc;  logg.modelres=[];  logg.modname=[];

% Do what?
logg.v1_1_ok={'p01_AF' ; 'p01_AS' ; 'p02_TY' ; 'p03_JV' ; 'p06_MK' ; 'p08_AK' ; 'p10_SG' ; 'p11_OG' ; 'p30_AR' ; 'p31_RD' ; 'p33_AA' ; 'p34_BE'};
logg.v1_2_ok={'p14_CP';'p15_SZ';'p23_SG';'p24_VP';'p26_PB';'p27_PV';'p29_MY';'p35_EO';'p37_ND';'p38_SL';'p39_JH';'p40_CC'; 'p45_SO';};
% eval(['logg.specificsubjects=logg.' logg.taskversion '_ok;'])
logg.taskversion='Both';  logg.specificsubjects=[];
logg.specificsubjects=[logg.v1_1_ok]; % logg.v1_2_ok];

% logg.specificsubjects=[logg.v1_1_ok];  logg.modelres='2b Hierar fits/res_hierarfitmodels_cF (16-Apr-2015) v1 bpjm10_L3066';               % v1 alone 
logg.specificsubjects=[logg.v1_2_ok]; logg.modelres='2b Hierar fits/res_hierarfitmodels_cF (16-Apr-2015) v2 bpj16_L3489';               % v1 alone 



% logg.specificsubjects= logg.specificsubjects([1:3 6:11]);


for o1=1:1 % General setup
    
    request.Data2Load= []; % ['(27-Feb-2014) Subjdata glmfit ' logg.taskversion]; % Empty to compile from scratch

    % Folders
    w=pwd; if strcmp(w(2), ':')==1;  where.modelfol='D:\Dropbox\SCRIPPS\8 Explore Mem\3b Modelling'; where.where='D:\Dropbox\SCRIPPS\8 Explore Mem'; else; 
        where.modelfol='/Users/EleanorL/Dropbox/SCRIPPS/8 explore mem/3b modelling'; 
        
        where.where='/Users/EleanorL/Dropbox/SCRIPPS/8 Explore Mem';  end
    where.data=[where.where filesep '2 All data' filesep logg.taskversion]; % Append version
    where.analysisscript=[where.where filesep '3a Analysis scripts'];
    where.analysisinputs=[where.where filesep '4 Inputs'];
    path(pathdef); addpath(where.analysisscript);  addpath(where.modelfol),  addpath([where.modelfol filesep '2 Val fxn']); addpath(genpath([where.modelfol filesep '2 Val fxn']))
    if isdir(where.analysisinputs)==0; mkdir(where.analysisinputs); end
    
    
    % Fetch subjects
    [w.s w.s1 logg.koshertable]=xlsread([where.analysisscript filesep 'i_subjectsok.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(logg.koshertable, logg.specificsubjects, logg.koshertable, 'Include');

    % Modelling results
    if isempty(logg.modelres)==0
        logg.mres=load([where.modelfol filesep '1 Inputs' filesep logg.modelres '.mat']);
        logg.modname=logg.mres.r_res{1,1};  logg.modpar=logg.mres.r_res{1,2};
        logg.mod=logg.mres.details.models(strcmp(logg.mres.details.models(:,1), logg.mres.r_res{1,1}),:);
        if isempty(strfind(logg.modname, 'w'))==0, logg.mod_explorepar=logg.modname(strfind(logg.modname, 'w')-1);  else, logg.mod_explorepar=[]; end
        fPar.cF_FL=  -8;
        fPar.cF_EC=  -2;
    else
        logg.mod_explorepar='no modelling requested';
    end
    
    disp(' ##################################################')
    disp(['No. of subjects:  ' num2str(logg.n_subjs)])
    if isempty(logg.specificsubjects)==0; disp(logg.subjects); end 
    disp(['Exploration-related variable:  ' logg.mod_explorepar])  % y= binomial variance, k=std dev
    input('Hit enter to continue     '); disp(' ')
    disp(' ##################################################')
end 

%% Fetch data (encoding and mem)
% subjdata: Col 1=Name, Col 2=Raw data & data in cells, Col 3= Data with just parameters needed for modelling ('col')

logg.rescale_variables=1;  % Re-scale IVs to [0 1]?
logg.DistortEnvThreat=0;
%
d_choice= [{'$' '$ explore' 'p(Accept|ExploreNoBomb)' 'p(Explore)'}; num2cell(nan(logg.n_subjs, 4))];  % Money won, Money won exploreonly, % Accept after explore-nobomb, % Explore
d_outcome=cell(logg.n_subjs+1,1); d_outcome{logg.n_subjs+1}=zeros(4,4);
if isempty(request.Data2Load)
    for o1=1:1 % Cols
        col.Subject=1;
        
        % Memory
        col.Roc=2;
        col.OldNew=3;
        col.RemKnow=4;
        col.SureGuess=5;
        col.Assoc=6;
        col.Mem=7; % Free floating memtype
        
        % Parameters
        col.Choice=8;
        col.EnvThreat=9;
        col.pLoss=10;
        col.NTokens=11;
        col.EV=12;
        col.Entropy=13;
        col.EntropyNTok=14;
        col.pLossSq=15;
        col.EntropyNTokSq=16;
        col.pLossOffset=17;
        col.EntropyNTokOffset=18;
        col.OutcomeMag=19;        
        col.OutcomeMean=20;
        col.OutcomeVariance=21;
        col.BinomVar=32;
        col.StanDev=33;
        col.EnvThreatOrig=34;
        col.VExplore=35;
        col.NTokensOrig=36;
        col.obs_pAccept=37;
        col.obs_pReject=38;
        col.obs_pExplore=39;
        col.vAccept=40;
        col.vReject=41;
        col.vExplore=42;
        col.Task=43;
        col.Trialnum=44;

        % Encoding performance
        col.encRT=30;
        col.encAcc=31;
        col.encRTmc=32;
            
        
    end
    subjdata=cell(logg.n_subjs+1,2);  subjdata{logg.n_subjs+1,1}='All'; subjdata{logg.n_subjs+1,3}.d=[]; subjdata{logg.n_subjs+1,3}.alld=[];
    for s=1:logg.n_subjs  % Fetch data from all subjects, and task parameters from 1st subject (col, design)
        subjdata{s,1}=logg.subjects{s};
        ws.enc=load([where.data filesep logg.subjects{s} filesep logg.subjects{s} '_file_3ExploreEnc.mat']);
        subjdata{s,2}.col.e=ws.enc.cF.col; subjdata{s,2}.encdata=ws.enc.cF.rdata; 
        subjdata{s,2}.design=ws.enc.cF.par.design; ws.encd=subjdata{s,2}.encdata; 
        ws.mem=load([where.data filesep logg.subjects{s} filesep logg.subjects{s} '_file_4Memtest.mat']);
        subjdata{s,2}.memdata=ws.mem.mem.data; subjdata{s,2}.col.m=ws.mem.mem.col;
        
        % Collate data (memdata)
        if s==1;  % columns
            col.e=ws.enc.cF.col;
            col.m=ws.mem.mem.col;
            col.m.EnvThreat_Level=21;
            col.m.NTokens_Level=22;
            col.m.NTokens=23;
            col.m.pLoss=24;
            col.m.Entropy=25;
            col.m.EntropyNTok=26;
            col.m.EV=27;
            col.m.OutcomeMean=28;
            col.m.OutcomeVariance=29;
            col.m.encRT=30;
            col.m.encAcc=31;
            %
            col.m.BinomVar=32;
            col.m.StanDev=33;
            col.m.EnvThreatOrig=34;
            col.m.VExplore=35; 
            
            % Insert new column/criteria here **
            
        end
        ws.new=subjdata{s,2}.memdata( subjdata{s,2}.memdata(:, col.m.Item_OldNew)==2,:) ;
        ws.new(:, [col.m.NTokens col.m.pLoss col.m.Entropy col.m.EntropyNTok col.m.EV col.m.OutcomeMean col.m.OutcomeVariance])=nan;
        ws.old=subjdata{s,2}.memdata( subjdata{s,2}.memdata(:, col.m.Item_OldNew)==1,:) ;
        ws.old(:,col.m.EnvThreat_Level)=ws.old(:,col.m.EnvThreat);
        for i=1:length( subjdata{s,2}.design.NTokenPairs_Npairs)
            ws.old( ws.old(:,col.m.NTokenPairs)==i, col.m.NTokens_Level)=i;
        end
        ws.old(:, col.m.NTokens)=ws.old(:, col.m.NTokens_Level).*2;
        for t=1:size(ws.old,1)  % Fetch encoding-stage variables
            wt.enc_trial=find(ws.old(t, col.m.ItemStim)==ws.encd(:,col.e.ItemStim));
            ws.old(t, col.m.encRT)=ws.encd(wt.enc_trial, col.e.RT1); 
            ws.old(t, col.m.encAcc)=ws.encd(wt.enc_trial, col.e.Trialvalid);
            wt=[];
        end
        
        
        % SORT
        ws.old= sortrows(ws.old, [col.EnvThreat col.NTokens]);
        
        % Write parameters into memdata
        [ ws.old ] = fpar_conflict( ws.old , col.m); 
        ws.old=sortrows(ws.old, [col.EnvThreat col.NTokens]);
        ws.new(:,[col.m.encRT col.m.encAcc])=nan;
        
        % Distort EnvThreat if needed + knock ons
        ws.old(:, col.m.EnvThreatOrig)= ws.old(:, col.m.EnvThreat);
        if logg.DistortEnvThreat==1 & isempty(strfind(logg.modname, 'j'))==0
            if s==1; input('Distort EnvThreat ok?'); end
            ss=find(strcmp(logg.mres.details.subjects,logg.subjects{s}));
            if isempty(strfind(logg.modname, 'f'))==0
                ws.mv.FixedLoss    = logg.modpar( ss,  3+find(strcmp(logg.mod{3}, 'f')));
            else  ws.mv=[];
            end
            ws.old(:, col.m.EnvThreat)=power(ws.old(:, col.m.EnvThreat),  logg.modpar( ss,  3+find(strcmp(logg.mod{3}, 'j')))); 
            [ ws.ov] = fcf_changeEnvThreat(ws.old(:, col.m.EnvThreat), ws.old(:, col.m.NTokens), ws.mv);
            ws.changevar={'EnvThreat';'NTokens';'pLoss';'Entropy';'EntropyNTok';'VExplore';'EV'};
            for c=1:length(ws.changevar)
                eval(['ws.old(:, col.m.' ws.changevar{c} ')=ws.ov.' ws.changevar{c} ';'])
            end
        end

        % Collate for modelling (translate columns from memdata space to modelling space)
        ws.mod(:,col.Roc)=ws.old(:, col.m.Roc);
        ws.mod(:,col.OldNew)=ws.old(:, col.m.OldNew);
        ws.mod(:,col.RemKnow)=ws.old(:, col.m.RemKnow);
        ws.mod(:,col.SureGuess)=ws.old(:, col.m.SureGuess);
        ws.mod(:,col.Assoc)=ws.old(:, col.m.CorrectPosition);
        ws.mod(:,col.EnvThreat)=ws.old(:, col.m.EnvThreat);
        ws.mod(:,col.EnvThreatOrig)=ws.old(:, col.m.EnvThreatOrig);
        ws.mod(:,col.pLoss)=ws.old(:, col.m.pLoss);
        ws.mod(:,col.NTokens)=ws.old(:, col.m.NTokens);
        ws.mod(:,col.NTokensOrig)=ws.mod(:,col.NTokens);
        ws.mod(:,col.EV)=ws.old(:, col.m.EV);
        ws.mod(:,col.Entropy)=ws.old(:, col.m.Entropy);
        ws.mod(:,col.Choice)=ws.old(:,col.m.Choice);
        ws.mod(:,col.EntropyNTok)=ws.old(:,col.m.EntropyNTok);
        ws.mod(:,col.encRT)=ws.old(:,col.m.encRT);
        ws.mod(:,col.encAcc)=ws.old(:,col.m.encAcc);
        ws.mod(:, col.OutcomeMag)=ws.old(:, col.m.OutcomeMag);
        ws.mod(:, col.BinomVar)=ws.old(:, col.m.BinomVar);
        ws.mod(:, col.StanDev)=ws.old(:, col.m.StanDev);
        
        
        % Get model values
        if isempty(logg.modelres)==0
            [ ws.transpar] = f_transpar(logg.mod{3}, logg.modpar(s, 4:end), 'from');
            ws.design(:,  [col.EnvThreat col.NTokens]) = ws.mod(:,  [col.EnvThreatOrig col.NTokensOrig]);
            ws.design(:,  col.Task    )=1;   ws.design(:,  col.Trialnum)=1:size(ws.design,1);
            eval(['[ ws.v] = ' logg.mod{1}  '(  ws.transpar,  {[] ws.design  fPar   col} ) ;'])
            ws.mod(:,  [col.vAccept  col.vReject  col.vExplore]) = squeeze(ws.v);
        else  ws.mod(:,  [col.vAccept  col.vReject  col.vExplore]) =nan;
        end
        
        
        
        for e=1:4 % Variables that belong to the 4x4 space 
            for n=1:4
                wc=ws.old(  ws.old(:, col.m.EnvThreatOrig)*4==e & ws.old(:, col.m.NTokens)/2==n, col.OutcomeMag);
                ws.mod( ws.mod(:, col.EnvThreatOrig)*4 ==e &  ws.mod(:, col.NTokens)/2 ==n,  col.OutcomeMean)=mean(wc);
                ws.mod( ws.mod(:, col.EnvThreatOrig)*4 ==e &  ws.mod(:, col.NTokens)/2 ==n,  col.OutcomeVariance)=var(wc);
                %
                wc=ws.mod( ws.mod(:, col.EnvThreatOrig)*4 ==e &  ws.mod(:, col.NTokens)/2 ==n,   col.Choice);
                ws.mod( ws.mod(:, col.EnvThreatOrig)*4 ==e &  ws.mod(:, col.NTokens)/2 ==n,  col.obs_pAccept)= mean(wc==1);
                ws.mod( ws.mod(:, col.EnvThreatOrig)*4 ==e &  ws.mod(:, col.NTokens)/2 ==n,  col.obs_pReject)= mean(wc==2);
                ws.mod( ws.mod(:, col.EnvThreatOrig)*4 ==e &  ws.mod(:, col.NTokens)/2 ==n,  col.obs_pExplore)= mean(wc==3);
            end
        end
        
        % New variables
        ws.mod(:,col.pLossSq)=ws.mod(:,col.pLoss).*ws.mod(:,col.pLoss);
        ws.mod(:,col.EntropyNTokSq)=ws.mod(:,col.EntropyNTok).*ws.mod(:,col.EntropyNTok);
        
        % Re-scale variables to [0 1]
        ws.mod_unscaled=ws.mod;
        if logg.rescale_variables
            logg.rescaled_vars={'EnvThreat';'NTokens';'pLoss';'Entropy';'EV'; 'EntropyNTok'; 
            'BinomVar';'StanDev';'OutcomeMean';'OutcomeMag';'OutcomeVariance';
            'encRT'; 'pLossSq';'EntropyNTokSq';
            'obs_pAccept';'obs_pReject';'obs_pExplore'
            };
            if s==1; disp('Rescaling variables to between 0 and 1'); disp(logg.rescaled_vars); end
            
            for v=1:length(logg.rescaled_vars)
                eval(['wv.col=col.' logg.rescaled_vars{v} ';'])
                wv.d=ws.mod(:,wv.col);
                
                if min(wv.d)<0, wv.d=wv.d+ abs(min(wv.d));   % Adjust floor
                elseif min(wv.d)>0, wv.d=wv.d - min(wv.d);
                end
                wv.d=wv.d./max(wv.d);            % Scale by ceiling
                %             disp([logg.rescaled_vars{v}  ':  '  num2str(min(wv.d)) ' to ' num2str(max(wv.d))])
                ws.mod(:,wv.col)=wv.d;
            end
        end
        
        % ################################################
        % Choice stats
        ws.e=ws.enc.cF.data;
        ws.eexp=ws.e(ws.e(:, col.e.Resp1)==3, :);
        ws.eexp_nb=ws.eexp(ws.eexp(:, col.e.BombExplored)==0, :);
        d_choice{s+1,1}=sum(ws.e(:, col.e.OutcomeMag));
        d_choice{s+1,2}=sum(ws.eexp(:, col.e.OutcomeMag));
        d_choice{s+1,3}=mean(ws.eexp_nb( :, col.e.Resp2)==1);
        d_choice{s+1,4}=mean(ws.e(:, col.e.Resp1)==3);
        
        %
        ws.mod(:,col.Subject)=s;
        subjdata{s,3}.alld=ws.mod;
        subjdata{s,3}.d=ws.mod(ws.mod(:,col.encAcc)==1, :);
        subjdata{s,3}.d_raw=ws.mod_unscaled(ws.mod_unscaled(:,col.encAcc)==1, :);  % with IVs unscaled
        subjdata{logg.n_subjs+1,3}.d=[subjdata{logg.n_subjs+1,3}.d; ws.mod];
        subjdata{logg.n_subjs+1,3}.alld=[subjdata{logg.n_subjs+1,3}.alld; subjdata{s,3}.alld];
        ws=[];
    end
    
    %
%     d_outcome{logg.n_subjs+1}=d_outcome{logg.n_subjs+1}./logg.n_subjs;
%     imagesc(d_outcome{logg.n_subjs+1}), colorbar, axis square
%     
%     imagesc(d_outcome{1}), colorbar, axis square
%     caxis([0.15 0.3])




    for o1=1:1 % Mark new variables (check columns don't overlap w old)
        col.Choice_Accept=21;
        col.Choice_Reject=22;
        col.Choice_Explore=23;
        
        disp(col) % Careful no overlap w old!!!!
        for s=1:logg.n_subjs+1
            ws.d=subjdata{s,3}.d;
            
            % Append new columns
            ws.d(:,col.Choice_Accept)=ws.d(:,col.Choice)==1;
            ws.d(:,col.Choice_Reject)=ws.d(:,col.Choice)==2;
            ws.d(:,col.Choice_Explore)=ws.d(:,col.Choice)==3;
            
            %
            subjdata{s,3}.d=ws.d;
            ws=[];
        end
    end
end

% Plot something?
plotsomething=0;
for o=1:1
    if plotsomething
        
        % Compile data
        d_outmean{logg.n_subjs+1}=zeros(4,4);
        for s=1:logg.n_subjs
            ws.d=subjdata{s,3}.alld;
            for e=1:4
                ee=5-e;
                for n=1:4
                    wc=ws.d( ws.d(:, col.EnvThreatOrig)*4 ==e &  ws.d(:, col.NTokensOrig)/2 ==n,:);
                    %
                    d_outmean{s}(ee,n)= mean(wc(:,  col.OutcomeMean) );
                end
            end
            d_outmean{logg.n_subjs+1}= d_outmean{logg.n_subjs+1} + d_outmean{s};
        end
        d_outmean{logg.n_subjs+1}= d_outmean{logg.n_subjs+1}./logg.n_subjs;
        
        % Plots
        plotwhat={'outmean'};
        f.plotcols=2;  f.figwidth= 1000; f.figheight=900;  f.fontsize=20; f.fontname='PT Sans Caption';
        f.subplot_VerHorz=[0.1 0.1]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.1 0.1]; k=1;
        figure('Name', 'Variables ', 'NumberTitle', 'off', 'Position', [100 200 f.figwidth f.figheight], 'Color', 'w');  k=1;
        for p=1:length(plotwhat)
            eval(['d_plot=d_' plotwhat{p} '{logg.n_subjs+1};']);
            
            % Plot
            subplot(ceil(length(plotwhat)/f.plotcols),f.plotcols, k)
            imagesc(d_plot); axis square; colorbar
            title(plotwhat{p} , 'FontSize', f.fontsize, 'FontName', f.fontname);
            set(gca, 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname);
            
        end
        
    end
end



%% Settings for GLMFIT

% Which IVs?
ivs={
'Subject';% 'Choice';
% 'Choice_Accept';
'Choice_Reject';'Choice_Explore';

'EnvThreat'; 


 'NTokens';
'pLoss';
'Entropy';  
'EV';  % Basics 


% [ Exploration-related variables ] ---------------------------------------
'EntropyNTok';
% 'BinomVar';  % y
% 'StanDev';
% 'ObspExplore';


% [ Outcome related ] ---------------------------------------
% 'OutcomeMag'; 
% 'OutcomeMean'; 
% 'OutcomeVariance'
};

% [ Values, probaility] ###############################################

% ivs={'Subject'; 'Choice'
% %     'obs_pAccept';'obs_pReject';'obs_pExplore';
%     'vAccept'; 'vReject'; 'vExplore'
%     
% %     'Choice_Accept';'Choice_Reject';'Choice_Explore';
%     };




% Set up IVs
logg.iv_cols=[]; for i=1:length(ivs); logg.iv_cols=[logg.iv_cols 'col.' ivs{i} ' ']; end; logg.iv_cols=eval(['[' logg.iv_cols ']']); 
logg.ivs_indiv=ivs(2:end); logg.iv_cols_indiv=logg.iv_cols(2:end); % For individual subjects


% Memory types
logg.memtypes={'hit';'surehit';'rem';'know';'surerem';'sureknow';'assoc';}; 
logg.memtypes_names={'Recognition memory (hit rate)';
    'Confident recognition hits';'Remember';'Know';'Sure remember';'Sure know';'Position recall accuracy'}; % ;'roc'};
% 
% logg.memtypes={'hit';'assoc';'rem';'know'}; 
% logg.memtypes_names={'Recognition memory (hit)';'Position recall accuracy';'Remember';'Know'}; 


%% Mixed effects (1st level, 2nd leve) 

% Data variables: row=sub, col= Factors (see ivs for list)
%         dc_ = coefficient, dp_ = p vals, dsc_ = sig coefficients
for m=1:length(logg.memtypes);
    eval(['dc_' logg.memtypes{m} '= nan(logg.n_subjs, length(ivs));'])
    eval(['dsc_' logg.memtypes{m} '= nan(logg.n_subjs, length(ivs));'])
    eval(['dp_' logg.memtypes{m} '= nan(logg.n_subjs, length(ivs));'])
end
logg.pthreshold_indiv=0.1; 
dc_eRT= nan(logg.n_subjs, length(ivs)); dsc_eRT= nan(logg.n_subjs, length(ivs)); dp_eRT= nan(logg.n_subjs, length(ivs)); dc_eAcc= nan(logg.n_subjs, length(ivs)); dsc_eAcc= nan(logg.n_subjs, length(ivs));dp_eAcc= nan(logg.n_subjs, length(ivs));


% Fit first-level GLM for all subjects separately
for s=1:logg.n_subjs
    ws.d=subjdata{s,3}.d; 
    ws.alld=subjdata{s,3}.alld; 
    
    for m=1:length(logg.memtypes)
        disp([logg.subjects{s} '  -  ' logg.memtypes{m}])
        % Which memory type?
        %   Roc: 1=Sure New
        %           2= Unsure New
        %           3= U K
        %           4= S K
        %           5= U R
        %           6= S R
        switch logg.memtypes{m}
            case 'hit'; wm.d= ws.d(:,col.Roc)>2.5;   % Hit
            case 'surehit'; wm.d= ws.d(:,col.Roc)==4 | ws.d(:,col.Roc)==6 ;  % SureHit
            case 'rem'; wm.d= ws.d(:,col.Roc)==5 | ws.d(:,col.Roc)==6 ;  % Rem
            case 'know'; wm.d= ws.d(:,col.Roc)==3 | ws.d(:,col.Roc)==4 ;  % Know
            case 'surerem'; wm.d= ws.d(:,col.Roc)==6 ; % SureRem
            case 'sureknow'; wm.d= ws.d(:,col.Roc)==4 ; % SureKnow
            case 'assoc'; wm.d= ws.d(:,col.Assoc); % Assoc
            case 'roc'; wm.d= ws.d(:,col.Roc)==4;
            otherwise; error(['Unrecognized memory type: ' logg.memtypes{m}]);
        end
        
        % Fit memory!
        [b d stats]=glmfit(ws.d(:,logg.iv_cols_indiv), wm.d, 'binomial');
        ws.pFit= glmval(b, ws.d(:,logg.iv_cols_indiv), 'logit');
        d_ll(s,m)= sum(log(binopdf(wm.d, size(wm.d,1), ws.pFit)));
        
        if strcmp(logg.memtypes{m}, 'roc')==1
            [b d stats]=glmfit(ws.d(:,logg.iv_cols_indiv), wm.d, 'normal');
        end
        wm.dc=b'; wm.dp=stats.p';
        wm.dsc=wm.dc; wm.dsc(wm.dp>logg.pthreshold_indiv)=nan;
        %
        
        
        
        % Log
        eval(['dc_' logg.memtypes{m} '(s,:)=wm.dc;'])
        eval(['dp_' logg.memtypes{m} '(s,:)=wm.dp;'])
        eval(['dsc_' logg.memtypes{m} '(s,:)=wm.dsc;'])
        wm=[]; ws.d(:,col.Mem)=nan;
    end
%     
%     % Fit encoding-stage performance
%     [b d stats]=glmfit(ws.d(:, logg.iv_cols_indiv), ws.d(:, col.encRT)); % RT
%     dc_eRT(s,:)=b';dp_eRT(s,:)=stats.p'; 
%     ws.eb=b'; ws.ep=stats.p'; ws.eb(ws.ep>logg.pthreshold_indiv)=nan;
%     dsc_eRT(s,:)=ws.eb;
%     [b d stats]=glmfit(ws.alld(:, logg.iv_cols_indiv), ws.alld(:, col.encAcc)); % Accuracy
%     dc_eAcc(s,:)=b'; dp_eAcc(s,:)=stats.p';
%     ws.eb=b'; ws.ep=stats.p'; ws.eb(ws.ep>logg.pthreshold_indiv)=nan;
%     dsc_eAcc(s,:)=ws.eb;
    
    %
    ws=[];
end
r_2ndlevel=[[{' ' } logg.ivs_indiv'];   [logg.memtypes   cell(length(logg.memtypes),length(logg.ivs_indiv))]]; r_2ndlevel_stats=r_2ndlevel;
% Second level 
for m=1:length(logg.memtypes) 
    eval(['dc_mem=dc_' logg.memtypes{m} '(:, 2:end);']);
    [h p ci stats]=ttest(dc_mem);  % Stats % h=p<0.1; disp('Include trends!'); % Include trends?
    
%     % Mark significance
    r_2ndlevel(m+1,2:end)=num2cell(h);
    for i=1:length(h), if stats.tstat(i)<0 & h(i)==1; r_2ndlevel{m+1, 1+i}=-1; end; end
    
    % Report full stats
%     r_2ndlevel_stats{1,1}=['df = ' num2str(stats.df(1))];
    for hh=1:length(h) 
%         r_2ndlevel_stats{m+1, hh+1}=['t=' num2str(stats.tstat(hh)) ',  p='  num2str(p(hh),3) ',  mean='    num2str( mean(dc_mem(:, hh)),3 ) ',  sd=' num2str(stats.sd(hh),3)];
        
        if p(hh)<0.05
%             r_2ndlevel_stats{m+1, hh+1}=['t(' num2str(stats.df(1))  ')= ' num2str(stats.tstat(hh),3) ',  p='  num2str(p(hh),3)];
            r_2ndlevel_stats{m+1, hh+1}=['t= ' num2str(stats.tstat(hh),3) ',  p='  num2str(p(hh),3)];
        else r_2ndlevel_stats{m+1, hh+1}='-';
        end
    end
    
    
% %     r_2ndlevel(m+1,find(stats.tstat.*h<0)+1)=r_2ndlevel{m+1, find(stats.tstat.*h<0)+1}*-1;  %  Reflect direction
%     r_2ndlevel(m+1, find(p<0.1 & p>0.05)+1)=num2cell(repmat(0.5, 1, sum(p<0.1 & p>0.05)));
    
end
r_2ndlevel_stats =[r_2ndlevel_stats ;[{' ' } num2cell(1:size(r_2ndlevel_stats,2)-1)]];

openvar r_2ndlevel, openvar r_2ndlevel_stats
% error('done :)')
% mean(d_ll)  % Likelihoods?




% Plots
doplot2ndlev=1;
if doplot2ndlev
%     ivs_names={'Choice';'Env Threat';'N Tokens';'p(Loss)'; 'Uncertainty';'EV' ; 'Value-scaled Uncertainty'};
    %
    close all hidden
    f.plotcols=3;  f.figwidth= 1000; f.figheight=900;  f.fontsize=20; f.fontname='PT Sans Caption';
    f.subplot_VerHorz=[0.1 0.1]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.1 0.1]; k=1;
    figure('Name', 'Mean betas', 'NumberTitle', 'off', 'Position', [100 200 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for m=1:length(logg.memtypes)
        eval(['dc_mem=dc_' logg.memtypes{m} '(:, 2:end);']);
        
        % Plot
        subplot(ceil(length(logg.memtypes)/f.plotcols),f.plotcols, k)
        barwitherr(std(dc_mem)/sqrt(logg.n_subjs),  mean(dc_mem),'y'  ) ; % Mean +/- Std Error
        title(logg.memtypes_names{m}, 'FontSize', f.fontsize, 'FontName', f.fontname);
        %     set(gca, 'xtick', 1:length(logg.ivs_indiv),'xticklabel', logg.ivs_indiv, 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
        %     set(gca, 'xtick', 1:length(logg.ivs_indiv),'xticklabel', ivs_names, 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname)
        set(gca, 'FontSize', round(f.fontsize*0.8), 'FontName', f.fontname);
        %     ylim([-15 15])
        %     xlim([0    length(logg.ivs_indiv)+1])
        %     xticklabel_rotate
        ylabel('Mean parameter estimate')
        
        k=k+1;
    end
    
    % Legend
    w.offset=0.25;  f.fontsize_legend=15;
    subplot(ceil(length(logg.memtypes)/f.plotcols),f.plotcols, k); axis off 
    for t=1:length(logg.ivs_indiv)
        text(0, 1.3- (t-1)*w.offset, [num2str(t) ':  ' logg.ivs_indiv{t}], 'FontSize',f.fontsize_legend);
    end
    
end


error('Done second level')

%% BETA CORRELATIONS '''################################

% Compile betas for requested IVs + all memtypes
whichiv={'Entropy'; 'EntropyNTok'; 'EnvThreat'};
d_betas=[whichiv cell(length(whichiv),1)];    % d_betas{iv,2}(subject, memtype)
for iv=1:length(whichiv)
    wi.whichiv_no=find(strcmp(logg.ivs_indiv, whichiv{iv}));
    d_betas{iv,2}=nan(logg.n_subjs, length(logg.memtypes));
    for m=1:length(logg.memtypes)
        eval(['d_betas{iv,2}(:, m)=dc_' logg.memtypes{m} '(:, wi.whichiv_no );'])
    end
end


% Do mem effects correlated with modelling params?
whichmodpar='w';
corr_modmem={  %  (1) mod par, (2) mem-IV, mem-type (3) data (modpar, betas), (4) results
    whichmodpar  {'EntropyNTok' 'hit'};           
    whichmodpar  {'EntropyNTok' 'surehit'};
    whichmodpar  {'EntropyNTok' 'rem'};
    whichmodpar  {'EntropyNTok' 'surerem'};
    whichmodpar  {'EntropyNTok' 'assoc'};
	};            
for c=1:size(corr_modmem,1)
    wc.allmembetas=d_betas{strcmp(d_betas(:,1), corr_modmem{c,2}{1}),2};
    wc.memnum=strcmp(logg.memtypes, corr_modmem{c,2}{2});
    %
    wc.membetas= wc.allmembetas(:, wc.memnum);
    wc.modpars=logg.modpar(:, 3+ find(strcmp(logg.mod{3}, corr_modmem{c,1})));
    %
    corr_modmem{c,3}=[wc.modpars wc.membetas];
    
    [r p]= corr(wc.membetas, wc.modpars);
%     [r p]= corr(wc.membetas, wc.modpars, 'type', 'Kendall'); disp('Kendalls rho!'); 
    wc.stat=r; 
    
    if p<0.05, corr_modmem{c,4} = wc.stat;
        disp([ num2str(c) '  ' corr_modmem{c,1} ' & '  corr_modmem{c,2}{1} '-fx on '  corr_modmem{c,2}{2} ' :   r='  num2str(r,2) '     p= ' num2str(p)]);
    elseif p<0.1, corr_modmem{c,4} = wc.stat;
        disp([ num2str(c) '  ' corr_modmem{c,1} ' & '  corr_modmem{c,2}{1} '-fx on '  corr_modmem{c,2}{2} ' :   r='  num2str(r,2) '     p= ' num2str(p)]);
    end
end

doplot=0;
if doplot
    
    cornum=3;
    
    
    close all hidden; 
    figure('color','w');   f.FontSize=30; f.FontName='PT Sans';
    
    
    wc.d= corr_modmem{cornum, 3};
    wc.memname= logg.memtypes_names{ strcmp(logg.memtypes, corr_modmem{cornum,2}{2})};
    
    
    % DISTORTIONS? 
%     wc.d(:,1)=  wc.d(:,1)- min(wc.d(:,1));
    
    % wc.d=sortrows(wc.d,1);
    % wc.d=wc.d(1:11,:);
    % wc.d=sortrows(wc.d,2);
    % wc.d=wc.d(1:10,:);
    %
    % wc.d(:,1)= zscore(wc.d(:,1));
    % wc.d(:,2)= zscore(wc.d(:,2));
    
    
    
    [r p]=corr(wc.d); r, p
    % [r p]=corr(wc.d, 'type', 'Kendall'); r, p
    scatter(wc.d(:,1), wc.d(:,2) ), l=lsline; set(l, 'LineWidth',4)
    % xlabel('Model parameter (z-scored)', 'FontSize', f.FontSize, 'FontName', f.FontName);
    % ylabel(['Parameter estimate for ' wc.memname ' (z-scored)'],'FontSize', f.FontSize, 'FontName', f.FontName);
    xlabel('w parameter', 'FontSize', f.FontSize, 'FontName', f.FontName);
    ylabel(sprintf(['Parameter estimates \n for ' wc.memname 'memory' ]),'FontSize', f.FontSize, 'FontName', f.FontName);
    set(gca, 'FontSize', f.FontSize, 'FontName', f.FontName)
    
end




% Do u and o effects tradeoff against each other?
disp('####### Are u and o betas trading off? ######')
for m=1:length(logg.memtypes)
    [r p]=corr(d_betas{strcmp(d_betas(:, 1), 'Entropy'),2}(:, m),d_betas{strcmp(d_betas(:, 1), 'EntropyNTok'),2}(:, m));
    if p<0.1; disp([logg.memtypes{m} ':  r=' num2str(r) ',    p= ' num2str(p,3)])
    else disp([logg.memtypes{m} ':   nsf '])
    end
end



% Do IV effects correlate with choice stats
r_betachoice=[whichiv repmat( {[ [{' '}; logg.memtypes] [d_choice(1,:); cell(length(logg.memtypes), size(d_choice,2))]]} ,length(whichiv),1)];
for iv=1:length(whichiv)
    disp([whichiv{iv} '  ----- ' ])
    for m=1:length(logg.memtypes)
        for c=1:size(d_choice,2)
            [r p]=corr(d_betas{iv,2}(:, m),  cell2mat(d_choice(2:end, c)));
            
%             [r p]=corr(d_betas{iv,2}(:, m),  cell2mat(d_choice(2:end, c)),'type', 'Kendall');  
%             disp('Kendalss')
            
            
            if p<.05; 
                disp([logg.memtypes{m} '  &  ' d_choice{1,c} '  :  r='  num2str(r,2) ',   p='  num2str(p) ' * '])
                r_betachoice{iv,2}{m+1,c+1}=r;
            else  % disp([logg.memtypes{m} '   &  ' d_choice{1,c} '  :   nsf'])
                disp([logg.memtypes{m} '  &  ' d_choice{1,c} '  :  r='  num2str(r,2) ',   p='  num2str(p)])
%                 r_betachoice{iv,2}{m+1,c+1}=p;
            end
        end
    end
end
openvar r_betachoice





%% How accurate are these predictors at the first level?


dothis=0;
if dothis
    r_predictmem=nan(logg.n_subjs,  length(logg.memtypes));   % Subjects, memtype
    col.phit=30;
    col.psurehit=col.phit+1;
    col.prem=col.phit+2;
    col.pknow=col.phit+3;
    col.psurerem=col.phit+4;
    col.psureknow=col.phit+5;
    col.passoc=col.phit+6;
    for s=1:logg.n_subjs
        ws.dd=subjdata{s,3}.d;
        ws.dd_raw=subjdata{s,3}.d_raw;  % Get unscaled IVs, for sorting only
        
        % Separate into cells
        ws.d=cell(4,4);
        %     ws.dsumm=nan(4*4,  size(ws.dd,2)); k=1;
        for e=1:4
            for n=1:4
                wc.d=ws.dd(ws.dd_raw(:, col.EnvThreat)*4==e & ws.dd_raw(:, col.NTokens)/2==n, :);
                wc.d(:, col.phit)=mean(wc.d(:, col.OldNew)==1);
                wc.d(:, col.psurehit)=mean(wc.d(:, col.OldNew)==1 & wc.d(:, col.SureGuess)==1);
                wc.d(:, col.prem)=mean(wc.d(:, col.RemKnow)==1);
                wc.d(:, col.psurerem)=mean(wc.d(:, col.RemKnow)==1 & wc.d(:, col.SureGuess)==1);
                wc.d(:, col.pknow)=mean(wc.d(:, col.RemKnow)==2);
                wc.d(:, col.psureknow)=mean(wc.d(:, col.RemKnow)==2 & wc.d(:, col.SureGuess)==1);
                wc.d(:, col.passoc)=mean(wc.d(:, col.Assoc)==1);
                
                
                
                STOPPED HERE. whats a good way to do this?
                
                ws.d{e,n}=wc.d;
            end
        end
        ws.ivs=ws.dd(:, logg.iv_cols_indiv);
        for m=1:length(logg.memtypes)
            switch logg.memtypes{m} % Load observed memory
                case 'hit'; wm.d= ws.dd(:,col.Roc)>2.5;   % Hit
                case 'surehit'; wm.d= ws.dd(:,col.Roc)==4 | ws.dd(:,col.Roc)==6 ;  % SureHit
                case 'rem'; wm.d= ws.dd(:,col.Roc)==5 | ws.dd(:,col.Roc)==6 ;  % Rem
                case 'know'; wm.d= ws.dd(:,col.Roc)==3 | ws.dd(:,col.Roc)==4 ;  % Know
                case 'surerem'; wm.d= ws.dd(:,col.Roc)==6 ; % SureRem
                case 'sureknow'; wm.d= ws.dd(:,col.Roc)==4 ; % SureKnow
                case 'assoc'; wm.d= ws.dd(:,col.Assoc); % Assoc
                case 'roc'; wm.d= ws.dd(:,col.Roc)==4;
                otherwise; error(['Unrecognized memory type: ' logg.memtypes{m}]);
            end
            
            % Load betas for resquested memory type
            eval(['dc_mem=dc_' logg.memtypes{m} ';']);
            wm.subbetas=dc_mem(s, :);
            
            
            
            %         error
            %
            %
            STOPPED HERE. Whats the right way to do this?
            %
            %
            %         glmval(wm.subbetas',  ws.ivs, 'logit')
            %
            %         dc_mem(s, 1)
            %
            %         dc_mem(s, 2:end)
            %
            %         ws.ivs
            %
            %
            %
            %         wm.d
            %         error
            %
            %         ws.ivs
            
            
            
            
            
            
            col
            
        end
    end
end

%% Fixed effects model (betas in variable pres, pvalues in pvals)
% Encoding-related variables: acc, accp, rt, rtp

error(':) :) DONE with all analysis. Stopped at where abigails''s analysis starts');

fixedfx=0;
if fixedfx
    s=length(logg.specificsubjects)+1; ws.d=subjdata{s,3}.d; % Group data
    pres=cell(length(logg.iv_cols)+1, length(logg.memtypes)+1); pres(1,2:end)=logg.memtypes; pres(2:end,1)=ivs;pvals=pres;  kon=pres(1:2, :); kon{2,1}='Constant'; bs=pres;
    for m=1:length(logg.memtypes)
        switch logg.memtypes{m}% What memory criteria (DVs)?
            case 'hit'; ws.d(:,col.Mem)= ws.d(:,col.Roc)>2.5;   % Hit
            case 'surehit'; ws.d(:,col.Mem)= ws.d(:,col.Roc)==4 | ws.d(:,col.Roc)==6 ;  % SureHit
            case 'rem'; ws.d(:,col.Mem)= ws.d(:,col.Roc)==5 | ws.d(:,col.Roc)==6 ;  % Rem
            case 'know'; ws.d(:,col.Mem)= ws.d(:,col.Roc)==3 | ws.d(:,col.Roc)==4 ;  % Know
            case 'surerem'; ws.d(:,col.Mem)= ws.d(:,col.Roc)==6 ;  % SureRem
            case 'sureknow'; ws.d(:,col.Mem)= ws.d(:,col.Roc)==4 ;  % SureKnow
            case 'assoc'; ws.d(:,col.Mem)=ws.d(:,col.Assoc);    % Assoc
            case 'roc'; ws.d(:,col.Mem)= ws.d(:,col.Roc);
        end
        
        % Fit memory!
        [b d stats]=glmfit(ws.d(:, logg.iv_cols), ws.d(:, col.Mem),'binomial');
        res=[[{'constant'}; ivs]  num2cell(stats.p) num2cell(b)];
        openvar res;  disp([logg.memtypes{m} ' #############################']); disp(res); disp(' ');
        %     input('Continue?   '); disp(' ')
        
        %    % Write to results
        kon{2, m+1}=res{1,3};
        %    kon{i, m+1}=res{1,3}
        for i=2:size(res,1)
            if res{i,2}<0.01;
                pres{i, m+1}=[num2str(res{i,3},3) ' **'];
                pvals{i, m+1}=num2str(res{i,2},3);
                bs{i, m+1}=res{i,3};
            elseif res{i,2}<0.05;
                pres{i, m+1}=[num2str(res{i,3},3) ' *'];
                pvals{i, m+1}=num2str(res{i,2},3);
                bs{i, m+1}=res{i,3};
            elseif res{i,2}<0.1;
                pres{i, m+1}=num2str(res{i,3}, 2);
                pvals{i, m+1}=num2str(res{i,2},3);
                bs{i, m+1}=res{i,3};
            else
                pvals{i, m+1}=[];
            end
            
        end
        
    end
%     openvar pres; openvar pvals; error('See results, fit GLM done :)')
    
    % Fit encoding-related variables (attentional effects)
    s=length(logg.specificsubjects)+1; ws.ad=subjdata{s,3}.alld; ws.d=subjdata{s,3}.d;
    acc=cell(length(logg.iv_cols)+1, 1); acc(2:end,1)=ivs; accp=acc; rt=acc; rtp=acc; m=1;
    [b d stats]=glmfit(ws.ad(:, logg.iv_cols), ws.ad(:, col.encAcc)); % accuracy
    eres=[[{'constant'}; ivs]  num2cell(stats.p) num2cell(b)]; % openvar eres
    for i=2:size(eres,1)
        
        if eres{i,2}<0.01;
            acc{i, m+1}=[num2str(eres{i,3},3) ' **'];
            accp{i, m+1}=num2str(eres{i,2},3);
        elseif eres{i,2}<0.05;
            acc{i, m+1}=[num2str(eres{i,3},3) ' *'];
            accp{i, m+1}=num2str(eres{i,2},3);
        elseif eres{i,2}<0.1;
            acc{i, m+1}=num2str(eres{i,3}, 2);
            accp{i, m+1}=num2str(eres{i,2},3);
        end
    end
    [b d stats]=glmfit(ws.d(:, logg.iv_cols), ws.d(:, col.encRTmc)); % rt
    eres=[[{'constant'}; ivs]  num2cell(stats.p) num2cell(b)];
    for i=2:size(eres,1)
        
        if eres{i,2}<0.01;
            rt{i, m+1}=[num2str(eres{i,3},3) ' **'];
            rtp{i, m+1}=num2str(eres{i,2},3);
        elseif eres{i,2}<0.05;
            rt{i, m+1}=[num2str(eres{i,3},3) ' *'];
            rtp{i, m+1}=num2str(eres{i,2},3);
        elseif eres{i,2}<0.1;
            rt{i, m+1}=num2str(eres{i,3}, 2);
            rtp{i, m+1}=num2str(eres{i,2},3);
        end
    end
    % openvar acc; openvar rt; openvar accp; openvar rtp
end
   

%% Correlate memfx and encoding fx 

% Correlation results: col = memtype, row = IV factor
w.corres=[[{' '}; logg.memtypes]  [ivs_indiv';  cell(length(logg.memtypes), length(ivs_indiv))]]';
rr_memrt=w.corres;  rp_memrt=w.corres;  rr_memacc=w.corres; rp_memacc=w.corres; 
rsr_memrt=w.corres; rsr_memacc=w.corres; % s= only print the significant correlation results

for m=1:length(logg.memtypes)
    eval(['wm.dc=dc_' logg.memtypes{m} ';'])

    for v=1:length(ivs_indiv)
        
        % Memory correlated w RT?
%         [wm.rr_memrt{v,1} wm.rp_memrt{v,1}]= corr(abs(wm.dc(:, v)), abs(dc_eRT(:, v)));
        [wm.rr_memrt{v,1} wm.rp_memrt{v,1}]= corr((wm.dc(:, v)), (dc_eRT(:, v)));
        if wm.rp_memrt{v,1}<0.1
            wm.rsr_memrt{v,1}=wm.rr_memrt{v,1};
        else
            wm.rsr_memrt{v,1}=[];
            wm.rp_memrt{v,1}=[];
        end
            
        % Memory correlated w enc Acc?
        [wm.rr_memacc{v,1} wm.rp_memacc{v,1}]= corr(abs(wm.dc(:, v)), abs(dc_eAcc(:, v)));
        if wm.rp_memacc{v,1}<0.1
            wm.rsr_memacc{v,1}=wm.rr_memacc{v,1};
        else
            wm.rsr_memacc{v,1}=[];
            wm.rp_memacc{v,1}=[];
        end
        
    end
    
    % Record results
    rr_memrt(2:end, m+1)=wm.rr_memrt; 
    rp_memrt(2:end, m+1)=wm.rp_memrt;
    
    
    rsr_memrt(2:end, m+1)=wm.rsr_memrt; 
    
    rr_memacc(2:end, m+1)=wm.rr_memacc; 
    rp_memacc(2:end, m+1)=wm.rp_memacc;
    rsr_memacc(2:end, m+1)=wm.rsr_memacc;
    
    wm=[];
end



%% Plot betas




% logg.memtypes={'hit';'surehit';'rem';'know';'surerem';'sureknow';'assoc'}; % ;'roc'};
%

par='pLossOffsetSq';
pnum=find(strcmp(ivs, par));
eval(['x=min(subjdata{1,3}.d(:,col.' par ')) :0.1: max(subjdata{1,3}.d(:,col.' par '));'])
plotcol={'r'; 'b';'g';[]; 'm';[]; 'c';'k';'y'};  % r=red, b= blue, g=green, m=magenta, c=cyan
close all hidden yy=nan(length(x), length(logg.memtypes));
for m=1:length(logg.memtypes)
    if isempty(bs{pnum+1, m+1})==0
        y=kon{2,1+m} +        bs{pnum+1, m+1}*x;
        plot(x,y, plotcol{m}, 'LineWidth',4 )
        yy(:,m)=y;
%     'LineWidth',2
%     'MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5])

        if m<length(logg.memtypes); hold on;end
        disp([logg.memtypes{m} '  ' plotcol{m}])
%         input('continue? ')
        %
    end
end
title(['Changes in memory scores as a function of ' par ' (range in task: ' num2str(min(x),1) ' to ' num2str(max(x),3) ')'])

%  change axes!
