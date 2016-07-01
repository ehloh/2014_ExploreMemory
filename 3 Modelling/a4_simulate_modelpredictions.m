% Plot model predictions (for comparison with data): Generate & plot choice contingencies PREDICTED by (specified) models
clear all; close all hidden; clc
% clear all; clc

% Request
request.predchoiceplot_individuals=0; 
request.predchoiceplot_groupmeans=1;


for o=1:1 % Request models 
    details.modelfams.b={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow';'b17_yw';'b18_fyw';'b19_eyw';'b20_feyw'; 'b21_kw';'b22_fkw';'b23_ekw';'b24_fekw'};
    details.modelfams.bp={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow';'bp17_yw';'bp18_fyw';'bp19_eyw';'bp20_feyw'; 'bp21_kw';'bp22_fkw';'bp23_ekw';'bp24_fekw'};
    details.modelfams.bm={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow';'bm17_yw';'bm18_fyw';'bm19_eyw';'bm20_feyw'; 'bm21_kw';'bm22_fkw';'bm23_ekw';'bm24_fekw'};
    details.modelfams.bpm ={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow';'bpm17_yw';'bpm18_fyw';'bpm19_eyw';'bpm20_feyw'; 'bpm21_kw';'bpm22_fkw';'bpm23_ekw';'bpm24_fekw'};
    %
    details.modelfams.bi={'bi01';'bi02_f';'bi03_e';'bi04_uw';'bi05_vw';'bi06_ow';'bi07_fe';'bi08_fuw';'bi09_fvw';'bi10_fow';'bi11_euw';'bi12_evw';'bi13_eow';'bi14_feuw';'bi15_fevw';'bi16_feow';'bi17_yw';'bi18_fyw';'bi19_eyw';'bi20_feyw'; 'bi21_kw';'bi22_fkw';'bi23_ekw';'bi24_fekw'};
    details.modelfams.bpi ={'bpi01';'bpi02_f';'bpi03_e';'bpi04_uw';'bpi05_vw';'bpi06_ow';'bpi07_fe';'bpi08_fuw';'bpi09_fvw';'bpi10_fow';'bpi11_euw';'bpi12_evw';'bpi13_eow';'bpi14_feuw';'bpi15_fevw';'bpi16_feow';'bpi17_yw';'bpi18_fyw';'bpi19_eyw';'bpi20_feyw'; 'bpi21_kw';'bpi22_fkw';'bpi23_ekw';'bpi24_fekw'};
    details.modelfams.bmi ={'bmi01';'bmi02_f';'bmi03_e';'bmi04_uw';'bmi05_vw';'bmi06_ow';'bmi07_fe';'bmi08_fuw';'bmi09_fvw';'bmi10_fow';'bmi11_euw';'bmi12_evw';'bmi13_eow';'bmi14_feuw';'bmi15_fevw';'bmi16_feow';'bmi17_yw';'bmi18_fyw';'bmi19_eyw';'bmi20_feyw'; 'bmi21_kw';'bmi22_fkw';'bmi23_ekw';'bmi24_fekw'};
    details.modelfams.bpmi ={'bpmi01';'bpmi02_f';'bpmi03_e';'bpmi04_uw';'bpmi05_vw';'bpmi06_ow';'bpmi07_fe';'bpmi08_fuw';'bpmi09_fvw';'bpmi10_fow';'bpmi11_euw';'bpmi12_evw';'bpmi13_eow';'bpmi14_feuw';'bpmi15_fevw';'bpmi16_feow';'bpmi17_yw';'bpmi18_fyw';'bpmi19_eyw';'bpmi20_feyw'; 'bpmi21_kw';'bpmi22_fkw';'bpmi23_ekw';'bpmi24_fekw'};
    %
    details.modelfams.bj={'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow';'bj17_yw';'bj18_fyw';'bj19_eyw';'bj20_feyw'; 'bj21_kw';'bj22_fkw';'bj23_ekw';'bj24_fekw'};
    details.modelfams.bpj={'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow';'bpj17_yw';'bpj18_fyw';'bpj19_eyw';'bpj20_feyw'; 'bpj21_kw';'bpj22_fkw';'bpj23_ekw';'bpj24_fekw'};
    details.modelfams.bjm={'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow';'bjm17_yw';'bjm18_fyw';'bjm19_eyw';'bjm20_feyw'; 'bjm21_kw';'bjm22_fkw';'bjm23_ekw';'bjm24_fekw'};
    details.modelfams.bpjm={'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm09_fvw';'bpjm10_fow';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow';'bpjm17_yw';'bpjm18_fyw';'bpjm19_eyw';'bpjm20_feyw'; 'bpjm21_kw';'bpjm22_fkw';'bpjm23_ekw';'bpjm24_fekw'};
    %
    details.modelfams.bji={'bji01';'bji02_f';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji07_fe';'bji08_fuw';'bji09_fvw';'bji10_fow';'bji11_euw';'bji12_evw';'bji13_eow';'bji14_feuw';'bji15_fevw';'bji16_feow';'bji17_yw';'bji18_fyw';'bji19_eyw';'bji20_feyw'; 'bji21_kw';'bji22_fkw';'bji23_ekw';'bji24_fekw'};
    details.modelfams.bpji={'bpji01';'bpji02_f';'bpji03_e';'bpji04_uw';'bpji05_vw';'bpji06_ow';'bpji07_fe';'bpji08_fuw';'bpji09_fvw';'bpji10_fow';'bpji11_euw';'bpji12_evw';'bpji13_eow';'bpji14_feuw';'bpji15_fevw';'bpji16_feow';'bpji17_yw';'bpji18_fyw';'bpji19_eyw';'bpji20_feyw'; 'bpji21_kw';'bpji22_fkw';'bpji23_ekw';'bpji24_fekw'};
    details.modelfams.bjmi={'bjmi01';'bjmi02_f';'bjmi03_e';'bjmi04_uw';'bjmi05_vw';'bjmi06_ow';'bjmi07_fe';'bjmi08_fuw';'bjmi09_fvw';'bjmi10_fow';'bjmi11_euw';'bjmi12_evw';'bjmi13_eow';'bjmi14_feuw';'bjmi15_fevw';'bjmi16_feow';'bjmi17_yw';'bjmi18_fyw';'bjmi19_eyw';'bjmi20_feyw'; 'bjmi21_kw';'bjmi22_fkw';'bjmi23_ekw';'bjmi24_fekw'};
    details.modelfams.bpjmi={'bpjmi01';'bpjmi02_f';'bpjmi03_e';'bpjmi04_uw';'bpjmi05_vw';'bpjmi06_ow';'bpjmi07_fe';'bpjmi08_fuw';'bpjmi09_fvw';'bpjmi10_fow';'bpjmi11_euw';'bpjmi12_evw';'bpjmi13_eow';'bpjmi14_feuw';'bpjmi15_fevw';'bpjmi16_feow';'bpjmi17_yw';'bpjmi18_fyw';'bpjmi19_eyw';'bpjmi20_feyw'; 'bpjmi21_kw';'bpjmi22_fkw';'bpjmi23_ekw';'bpjmi24_fekw'};
    %
    details.ct_modnums=[1 3 4 5 6 11 12 13 17 19 21 23]; % true models for ct
end

% Hierarchical fits ##############
request.fittype='hierarfit';    request.results_folder='2b Hierar fits/';  
details.fitresults_date='(16-Apr-2015) v1 bpjm10_L3066'; details.dataset_date='v1 (01-Apr-2015)'; details.whichmodels=details.modelfams.bpjm(10);  
% details.fitresults_date='(16-Apr-2015) v2 bpj16_L3489';  details.dataset_date='v2 (01-Apr-2015)'; details.whichmodels=details.modelfams.bpj(16);  





% Fixed fits ##############
% request.fittype='fit';    request.results_folder='2a Fixed fits/';  
% details.fitresults_date='(05-Apr-2015) v1 bpjm18y_L3029';   details.dataset_date='v1 (01-Apr-2015)';
% details.fitresults_date='(05-Apr-2015) v2 bpjm16_L3428';   details.dataset_date='v2 (01-Apr-2015)';
% details.whichmodels=details.modelfams.bpji(16);  

for o1=1:1 %% Set up
    
    
    request.tasktype=1; % 1=cF, 2=ct
    % Folders
    w.where=pwd;
    if strcmp(w.where(2), ':')==1;   where.analysisscripts='D:\Dropbox\SANDISK\8 Explore Mem\3a Analysis scripts'; where.scripts='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling';
    else where.analysisscripts='/Users/EleanorL/Dropbox/sandisk/8 Explore Mem/3a Analysis scripts'; where.scripts='/Users/EleanorL/Dropbox/sandisk/8 Explore Mem/3b Modelling';        
    end
    cd(where.scripts); path(pathdef);   addpath(where.analysisscripts)
    request.tasktype_duringfit=request.tasktype; %  1=cF, 2=ct, 3=cFct
    details.modelmodfits=['res_' request.fittype 'models_cF ' details.fitresults_date  '.mat'];
%     task='conflict';  
    
    w.dd=load([where.scripts filesep '1 Inputs' filesep 'All data ' details.dataset_date '.mat']); subjdata=w.dd.subjdata; 
    d_fits=load([where.scripts filesep '1 Inputs' filesep request.results_folder details.modelmodfits]);
    details.subjects=subjdata(:,1); details.n_subjs=size(subjdata,1);
    col=rmfield(d_fits.details.col, {'OutcomeMagnitude'});
    details.fixedpar=d_fits.details.fixedpar;
    addpath([where.scripts filesep '2 Val fxn']); addpath(genpath('2 Val fxn'))
    
    % Model settings + fetch requested model details (in requested order)
    d_fits.details4fit=d_fits.details; d_fits=rmfield(d_fits, 'details');
    orig.r_res=d_fits.r_res; 
    orig.r_iterations=[]; % disp('r_iterations turned off! all other lines w r_iterations also turned off. search r_iterations')
    orig.r_iterations=d_fits.r_iterations; % change to grid if needed
    if isempty(details.whichmodels); details.whichmodels=orig.r_res(1:5,1); end
    details.n_models=length(details.whichmodels); 
    d_fits.r_res=cell(length(details.whichmodels),  size(d_fits.r_res,2));
    d_fits.r_iterations=cell(length(details.whichmodels),  size(d_fits.r_iterations,2));
    for m=1:  length(details.whichmodels)
        if sum(strcmp(orig.r_res(:,1), details.whichmodels{m}))==1
            wm.rownum=find(strcmp(orig.r_res(:,1), details.whichmodels{m}));
            d_fits.r_res(m, 1:length(orig.r_res(wm.rownum,:)))= orig.r_res(wm.rownum,:);
            %
            wm.rownum=find(strcmp(orig.r_iterations(:,1), details.whichmodels{m}));
            d_fits.r_iterations(m, 1:length(orig.r_iterations(wm.rownum,:)))= orig.r_iterations(wm.rownum,:);
            %
            wm=[];
        else error(['Cannot find requested model: '         details.whichmodels{m}])
        end
    end
    
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels);
    disp('====================================')
    disp(['Requested: ' num2str(size(subjdata,1)) ' subjects']); disp(' ')
    disp('Models requested:'); disp(details.whichmodels); disp(' ')
    % input('Hit enter to start                                   ');
    disp('====================================')
end

m=1;
details.subj_ntrials=cellfun(@(x) size(x,1),  subjdata(:,2)); L= sum(orig.r_res{m,2}(:,2));   R = sum(details.subj_ntrials) .* -log(1/3);
disp(['Pseudo r2 for winning model  - '  orig.r_res{m,1} '  :' num2str(   1-L/R )]);

% cellfun(@(x) strfind(x, 'j'), orig.r_res(:,1), 'UniformOutput',0)
% openvar ans

%% (1) Generated choices as predicted by models 
%   d_fits: results from model fit

% Data variables 
disp(['Winning model:    '  d_fits.r_res{1,1} ]);
d_pred=cell(details.n_subjs+1, details.n_models); % d_pred{sub,row}=simulated trialstats for each subject; 
p_choicematrix=cell(details.n_subjs+1, details.n_models,3); % p(Accept/Reject/Explore); row=sub, col=model, 3rd dim=choice
v_choicematrix=cell(details.n_subjs+1, details.n_models,3); % value(Accept/Reject/Explore); row=sub, col=model, 3rd dim=choice
d_bestchoice=cell(details.n_subjs+1, details.n_models,1);
d_fits.original_res=d_fits.r_res;

% Artificially alter details
% disp('Parameters altered!!')
% d_fits.r_res{1,2}=d_fits. r_res{1,2}(1,:); details.subjects=details.subjects(1); details.n_subjs=1; disp('One subject only!');
% d_fits.r_res{1,2}(1,3+1:end)=[1 0 1 -8 0 0];   % cF 
% p=5;  d_fits.r_res{1,2}(:,3+p)=0;
% d_fits.r_res{1,2}(:,4:end)=repmat([1 0 1   -12 0 0], 20,1);    % cF 



needtransform=0;
if needtransform
    input('Transform parameters?');
    mnum=1;
    modpars=d_fits.details4fit.models{strcmp(d_fits.details4fit.models(:,1), d_fits.r_res{mnum,1}),3};
    for s=1:size(d_fits.r_res{mnum,2},1)
        ws.pars= d_fits.r_res{mnum,2}(s, 4:end);
        
        % What transform?
        ws.newpars=f_transpar(modpars,  ws.pars, 'to');
        
        % Put back
        d_fits.r_res{mnum,2}(s, 4:end)=ws.newpars; ws=[];
    end
end

%% (1) Generate value & choices predicted by models


% Generate simulated choices from model predictions : d_pred{model, subject} & group
for o1=1:1
    
    % Columns (for data to be added)
    for o2=1:1
        col.evAccept=21;
        col.evReject=col.evAccept+1;
        col.evExplore=col.evAccept+2; % note: EV(Explore)=EV of choosing Explore; VExplore=EV(Explore) - EV(NonExplore)
        
        col.pAccept=24;
        col.pReject=col.pAccept+1;
        col.pExplore=col.pAccept+2;
        col.StanDev=col.pAccept+3;
        col.BinomVar=col.pAccept+5;
    end
    
    % Specify design variables
    d_design=nan(4*4, 10); d_design(:,  [col.EnvThreat col.NTokens])=[sortrows(repmat((1:4)',4,1))/4 2*repmat((1:4)',4,1)]; d_design(:,col.Task)=request.tasktype;
    [d_design] = fpar_conflict(d_design, col);
    
    for m= 1:details.n_models
        wm.model=details.whichmodels{m};
        wm.modrow_res=find(strcmp(d_fits.r_res, details.whichmodels{m}));
        %         disp(['Sum of nLLs: ' num2str(sum(d_fits.r_res{wm.modrow_res,2}(:,2)))])
        
        % Load mean parameters (group mean)
        for i=4:size(d_fits.r_res{wm.modrow_res,2},2)
            d_fits.r_res{wm.modrow_res,2}(details.n_subjs+1, i)=d_fits.r_res{wm.modrow_res, i+2};
        end
        
        % Subject-specific predictions followed by mean predictions
        for s=1:details.n_subjs+1
            ws.modpar=d_fits.r_res{wm.modrow_res,2}(s,4:end); % Fetch subject's model parameters
            ws.origpar=ws.modpar;
            
            % Apply inverse transforms
            for p=1:details.models{m,2}
                x=ws.modpar(p); eval(['ws.modpar(p)= ' details.models{m, 5}{p}  ';'])
            end
            
            
            if s<details.n_subjs+1
                if isempty(strfind(details.models{m,1}, 'p'))==1;
                    [nll pch]=f_nllsoftmax(ws.modpar, {details.models{m,1} subjdata{s, request.tasktype+1} details.fixedpar col});
                else
                    [nll pch]=f_nllsoftmax_lapse(ws.modpar, {details.models{m,1} subjdata{s, request.tasktype+1} details.fixedpar col});
                end
                if  d_fits.r_res{m,2}(s,2)-nll~=0;   disp(['nLL calc discrepancy!    '  details.models{m,1} '  (' subjdata{s,1} ')    '  num2str(d_fits.r_res{m,2}(s,2),4) ' vs ' num2str(nll, 4)]); end % Paranoid check
            end
            
            % V(Choices) according to this model
            d_pred{s,m}=d_design;
            eval(['[ d_pred{s,m}(:,[col.evAccept col.evReject col.evExplore])] = ' details.whichmodels{m} '(ws.modpar, {[] d_design details.fixedpar, col});'])
            
            % Calculate p(Choice) using softmax rule
            ws.beta=ws.origpar(1);
            ws.softmaxbase=( exp(ws.beta.*d_pred{s,m}(:,col.evAccept))+exp(ws.beta.*d_pred{s,m}(:,col.evReject))+exp(ws.beta.*d_pred{s,m}(:,col.evExplore))  );
            if isempty(strfind(details.models{m,1}, 'p'))~=1;
                ws.epsilon=ws.origpar(2);
                ws.softmaxbase= ws.softmaxbase .*(1-2.*ws.epsilon + ws.epsilon);
            end
            d_pred{s,m}(:,col.pAccept)=exp(ws.beta.*d_pred{s,m}(:,col.evAccept)) ./ ws.softmaxbase;
            d_pred{s,m}(:,col.pReject)=exp(ws.beta.*d_pred{s,m}(:,col.evReject)) ./ ws.softmaxbase;
            d_pred{s,m}(:,col.pExplore)=exp(ws.beta.*d_pred{s,m}(:,col.evExplore)) ./ ws.softmaxbase;
            
            % Format p(Choices) into 6x6 (for plotting)
            ws.d=d_pred{s,m};
            for e=1:4
                for n=1:4
                    
                    % pChoice
                    p_choicematrix{s, m, 1}(5-e,n)=ws.d(ws.d(:,col.EnvThreat) ==e/4 & ws.d(:,col.NTokens)==2*n, col.pAccept);
                    p_choicematrix{s, m, 2}(5-e,n)=ws.d(ws.d(:,col.EnvThreat) ==e/4 & ws.d(:,col.NTokens)==2*n, col.pReject);
                    p_choicematrix{s, m, 3}(5-e,n)=ws.d(ws.d(:,col.EnvThreat) ==e/4 & ws.d(:,col.NTokens)==2*n, col.pExplore);
                    
                    % V(Choice)
                    v_choicematrix{s, m, 1}(5-e,n)=ws.d(ws.d(:,col.EnvThreat) ==e/4 & ws.d(:,col.NTokens)==2*n, col.evAccept);
                    v_choicematrix{s, m, 2}(5-e,n)=ws.d(ws.d(:,col.EnvThreat) ==e/4 & ws.d(:,col.NTokens)==2*n, col.evReject);
                    v_choicematrix{s, m, 3}(5-e,n)=ws.d(ws.d(:,col.EnvThreat) ==e/4 & ws.d(:,col.NTokens)==2*n, col.evExplore);
                   
                    % Best choice
                    ws.vc=[v_choicematrix{s, m, 1}(5-e,n) v_choicematrix{s, m, 2}(5-e,n) v_choicematrix{s, m, 3}(5-e,n)];
                    if sum(ws.vc==max(ws.vc))==1
                        d_bestchoice{s, m,1}(5-e,n)=find(ws.vc==max(ws.vc));
                    elseif sum(find(ws.vc==max(ws.vc)))==4 % A + E
                        d_bestchoice{s, m,1}(5-e,n)=1.8;
                    elseif sum(find(ws.vc==max(ws.vc)))==5 % R + E
                        d_bestchoice{s, m,1}(5-e,n)=2.5;
                    elseif sum(find(ws.vc==max(ws.vc)))==3 % A + R
                        d_bestchoice{s, m,1}(5-e,n)=1.5;
                    end
                    
                end
            end
            
            % ######################################################################
            % INSERT HERE To use subject's model parameters to calculate other variables (theoretical best? Rank order choices?)
            
            
            
            
            
            
            % ######################################################################
            
            ws=[];
        end
        
        %
        wm=[];
    end
    
    % Re-calculate mean as the average of simulated individual data
    for m= 1:details.n_models
        wc.p=repmat({zeros(4,4)},3,1); wc.v=wc.p;
        for c=1:3
            
            % Compile
            for s=1:details.n_subjs
                wc.p{c}=wc.p{c}+  p_choicematrix{s, m, c};
                wc.v{c}=wc.v{c}+  v_choicematrix{s, m, c};
            end
            
            % Calculate mean
            p_choicematrix{details.n_subjs+1, m, c}= wc.p{c}/details.n_subjs;
            v_choicematrix{details.n_subjs+1, m, c}= wc.v{c}/details.n_subjs;
            
        end
    end
end





%% (1) Plot simulated values/choices for group & individuals

% plot WHAT?
plotmatrix=p_choicematrix; range=[0 1];
% plotmatrix=v_choicematrix; range=[-8 8];
% plotmatrix=d_bestchoice; range=[1 3];


% vv=plotmatrix(details.n_subjs+1,1:3);

close all hidden
% Plot predictions for overall group
if request.predchoiceplot_groupmeans==1;
% fontsize=15; fontname='PT Sans Caption';
    f.plotcols=4;   f.fontsize=15; f.fontsize_title=15;  f.fontname='PT Sans Caption';  
    f.subplot_VerHorz=[0.05 0.1]; f.fig_BotTop=[0.01 0.03]; f.fig_LeftRight=[0.02 0.02];
    figure('Name', 'Group means (Predicted choice)', 'NumberTitle', 'off', 'Position', [130 85 900 600], 'Color', 'w');  k=1;
    request.tasktype_choices={'% Accept'  '% NoBomb';  '% Reject'  '% Bomb';  '% Explore'  '% Explore'; };
    for m=1:details.n_models
        s=details.n_subjs+1; % Mean data is at the end of variable
        
        subtightplot(details.n_models,  f.plotcols, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        text(0.5,0.5 , details.whichmodels{m},'FontSize', f.fontsize, 'FontName', f.fontname);   axis 'off'; k=k+1;
        
        for c=1:size(plotmatrix,3)
            subtightplot(details.n_models,  f.plotcols, k, f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagesc(plotmatrix{s,m,c}); axis('square');colorbar;
            title(request.tasktype_choices{c,request.tasktype}, 'FontSize', f.fontsize_title, 'FontName', f.fontname)
            caxis(range); 
            
            
            % Labels
            ylabel('Environmental Threat','FontSize', f.fontsize, 'FontName', f.fontname);  
            xlabel('No. Tokens','FontSize', f.fontsize, 'FontName', f.fontname);  
            set(gca,'YTick',1:6, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick',1:4, 'XTickLabel',2:2:2*4)
            set(gca, 'FontSize', f.fontsize, 'FontName', f.fontname)
        end
    end
    
end







% Plot individuals' predicted choice
if request.predchoiceplot_individuals
    % Figure details
    f.figheight= 100; f.figheight=800; f.subplotcols=4; f.subplot_VerHorz=[0.01 0.07]; f.fig_BotTop=[0.01 0.01]; f.fig_LeftRight=[0.05 0.1];
    for m=1:details.n_models
        
        % Predicted choice
        fc=figure('Name', ['(Fig. ' num2str(m) ') ' details.whichmodels{m}  ' - Predicted individual choice (Accept, Reject, Explore)'], 'NumberTitle', 'off', 'Position',[200,00,f.figheight,f.figheight]); set(gcf,'Color',[1 1 1]);
        slist1=[1 3 5 7 9 11 13 15 20];
        slist2=[2 4 6 8 10 12 14 16 17 18 19];
        
%         details.n_subjs=length(slist1);
        for s=1:details.n_subjs
            % Label subjects
            subtightplot(details.n_subjs,f.subplotcols,(s-1)*f.subplotcols+1,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
            text(0.8,0.3, details.subjects{s}); axis 'off'
            for c=1:3 % For each choice
                subtightplot(details.n_subjs,f.subplotcols,(s-1)*f.subplotcols+1+c,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
                
                
                
                % Plot whatever is requested
                imagesc(plotmatrix{s,m,c});
                
                
                % Blah
                axis('square'); axis 'off' ; caxis(range); % colorbar
            end
        end
    end
end













