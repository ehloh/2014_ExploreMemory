% Fit models to the choice data, calculate BICs & parameter values
clear all; close all hidden; clc

% Requested analysis
details.maxgroupiter=200; % Max no. iters!
ThisPC='Mac';
% 
% details.dataset_date='both (01-Apr-2015)';
details.dataset_date='v1 (01-Apr-2015)';
% details.dataset_date='v1 (31-Mar-2015) original';
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
for o=1:1
    
    topatch.both={'bpjmi10_fow';'bpjmi11_euw';'bpjmi13_eow';'bpjmi14_feuw';'bpjmi15_fevw';'bpjmi16_feow';'bpjmi17_yw';'bpjmi18_fyw';'bpjmi19_eyw';'bpjmi20_feyw';'bpjmi21_kw';'bpjmi22_fkw';'bpjmi23_ekw';'bpjmi24_fekw';};
    topatch.v1={'bpjm24_fekw';'bpm10_fow';'bpmi17_yw';'bpmi18_fyw';'bpmi19_eyw';'bpmi20_feyw';'bpmi21_kw';'bpmi22_fkw';'bpmi23_ekw';'bpmi24_fekw';};
    topatch.v2={'bpjmi07_fe';'bpjmi08_fuw';'bpjmi09_fvw';'bpjmi10_fow';'bpjmi11_euw';'bpjmi12_evw';'bpjmi13_eow';'bpjmi14_feuw';'bpjmi15_fevw';};
    topatch.v1orig={'bi19_eyw';'bjmi19_eyw';'bjmi20_feyw';'bmi19_eyw';'bpi17_yw';'bpj05_vw';'bpj16_feow';'bpj24_fekw';'bpji17_yw';'bpji22_fkw';'bpji24_fekw';'bpjm07_fe';'bpjm14_feuw';'bpjmi15_fevw';'bpjmi17_yw';'bpmi10_fow';'bpmi22_fkw';};
    
    
    
end



which=1:24; % which =details.modelfams.ct_modnums;
details.whichmodels=[
%     details.modelfams.b(which)
%     details.modelfams.bp(which)
%     details.modelfams.bm(which)
%     details.modelfams.bpm(which)
% %     %
%     details.modelfams.bi(which)
%     details.modelfams.bpi(which)
%     details.modelfams.bmi(which)
%     details.modelfams.bpmi(which)
%     
    details.modelfams.bj(which)
    details.modelfams.bpj(which)
    details.modelfams.bjm(which)
    details.modelfams.bpjm(which)
%      %
%     details.modelfams.bji(which)
%     details.modelfams.bpji(which)
%     details.modelfams.bjmi(which)
%     details.modelfams.bpjmi(which)
    ];

% details.whichmodels=details.modelfams.bpjmi(14);

for o1=1:1 %% Set up
    
    details.tasktype=1; % 1=cF, 2=ct
    Valfxn_type=[];   details.folder4saving=[]; details.Valfxn_type=Valfxn_type; w=pwd;  
    if strcmp(w(2), ':')==1;   where.where='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling'; 
    else where.where='/Users/EleanorL/Dropbox/SANDISK/8 Explore Mem/3b Modelling'; 
    end
    path(pathdef);  addpath([where.where filesep '2 Val fxn'])
%     addpath(genpath('2 Val fxn'))
    addpath([where.where filesep '2 Val fxn' Valfxn_type]);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'b']);   % 
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bp']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bm']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpm'])
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bi']);   % 
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpi']); 
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bmi']); 
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpmi']);  
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bj']);     % 
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpj']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bjm']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpjm']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bji']);    % 
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpji']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bjmi']);
    addpath([where.where filesep '2 Val fxn' Valfxn_type filesep 'bpjmi']);
    if details.tasktype==1; details.task='cF'; elseif details.tasktype==2; details.task='ct'; else error('Invalid details.tasktype requested (1=cF, 2=ct)'); end    
    
    % Task parameters
    details.dataset=['All Data ' details.dataset_date '.mat']; % Selection of specific subjects not implemented
    w=load([where.where filesep '1 Inputs' filesep details.dataset]);
    col=w.details.col; subjdata=w.subjdata; details.n_subjs=w.details.n_subjs;  details.subjects=w.details.subjects;  details.col=col;
    details.fixedpar.cF_FL=  -8;
    details.fixedpar.cF_EC=  -2;
    details.fixedpar.ct_EC=  0;
        for s=1:details.n_subjs; details.subj_ntrials(s,1)=length(subjdata{s, details.tasktype+1});end;   % n trials per subject
        rand('state',sum(100*clock));
    details.resetZ=[];
    diary([where.where filesep '1 Inputs' filesep '1a Fit logs' filesep 'diary_fithierar_' details.task ' (' date ')'])
    
    % Fit/fminunc settings 
    details.tasktype=1; % 1=cF, 2=ct, 3=Both
    details.fitgroup_n_eviter=1;
    details.resetZ=[];
    details.fitsub_n_iter=2;
    details.calcbic_n_samples=10000;
    w.options=optimset('display','off','DerivativeCheck','on');
    warning('off','optim:fminunc:SwitchingMethod')
    
    % Model settings + fetch requwsi.fitparvalsed model details
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels);
    details.n_models=length(details.whichmodels);
    errorlog={}; e=1;
    details.modelgroupfitok=nan(details.n_models,1);
    
    % Model settings + fetch requested model details
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels,details.fitsub_n_iter);
    details.n_models=length(details.whichmodels); errorlog={}; e=1;
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
    
    disp('=============================================================='); w.c=clock; w.c1=w.c;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ');
    disp(['Requested: ' num2str(details.n_subjs) ' subjects'])
    disp(['details.task type:  ' details.task ]); disp(' ');
    disp([num2str(length(details.whichmodels))  ' Models:']); disp(details.whichmodels); disp(' ');
    disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
    disp(['Max no. iters = '  num2str(details.maxgroupiter)]); disp(' ')
    input('Hit enter to start                                   ');
    disp('==============================================================')
end

%% Run all models

details.do.dofit=1;
details.do.doreturn=0;


if details.do.dofit
    
    startfrom=1;
    we=[];
    
    % 'r_res'
    %       Col 1: Model name
    %       Col 2: Subject fit parameters
    %             Col i: BIC (omitted!)
    %             Col ii: nLL
    %             Col iii: fminunc exceeded default iterations?
    %             Col iv onwards: parameters (beta first)
    %       Col 3: Model BIC (summed across subjects)
    %       Col 4: Hessians
    %       Col 5:
    %       Col 6 onwards: mean parameter values
    for o1=1:1 % Columns
        rc.modname=1;
        rc.subpars=2;
        rc.sp.nll=2;
        rc.sp.p1=4;
        rc.modelbic=3;
        rc.hessians=4;
        rc.mean_p1=6;
        %
        
        if exist('r_res','var') ==0; r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2))));  r_iterations=cell(details.n_models,4); r_iterd=cell(details.n_models,2); end
    end
    for m=startfrom:details.n_models
        w.c=clock; disp(['Model ' num2str(m) ' - ' details.models{m,1} '       ['  num2str(w.c(4))  ':'  num2str(w.c(5))  '] ##############' ])
        
        % r_iterations -  model fits for each subject x iteration
        %       Col 1               = Model name
        %       Col 2 onwards  = Within-EV Iteration results
        %                                   [  1 row of r_iterations r{m,:}  =  {modelname  ev1_ParVariable   ev1_HessVariable ........  ev2_ParVariable   ..}  ]
        %                                           ev# =requested iteration of ev (outside all ev steps)
        %                 [ r_iterations{m,1+2*(i_evmod-1)+1}{within-ev-iter,1} ]
        %                         Col 1: Subject            [params]
        %                         Col 2: nLL actual
        %                         Col 3: posterior nLL (with population prior)
        %                         Col 4: N iterations
        %                         Col 5 onwards: model parameters (1st parameter is beta/inverse temperature)
        %                                               (Params are in fit-space right until separate correction at the end of all modelfits)
        %                 [ r_iterations{m,1+2*(i_evmod-1)+1}{within-ev-iter,2} ]
        %                         Col 1: Subject            [params]
        %                         Col 2: Hessians
        %                         Col 3: Variance
        %                                    (All values are in fit-space not param space)
        % r_iterd - details of each ev step fit (within each requested ev procedure)
        %       Col 1               Model name
        %       Col 2               Data for plotting convergence (row=subject, col=iteration)
        %       Col 3               Details for requested iter 1..
        r_iterations{m,1}=details.models{m,1};
        r_iterd{m,1}=details.models{m,1};
        r_res{m,1}=details.models{m,1};
        
        for i_evmod=1:details.fitgroup_n_eviter % Requested iterations of entire EV protocol
            disp(['Requested iteration #' num2str(i_evmod) ' --------------------------------------------------------'])
            
            %% Fit group data (find convergence  at group level, can take >1000 runs)
            %   i_group/wg= iteration within EV process (automated),  i_evmod/we=requested iteration of entire model fit (outside group EV protocol)
            i_group=1;  wg=[]; Z=[];
            while 1;
                r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}=[(1:details.n_subjs)'     nan(details.n_subjs, 3+details.models{m,2})];
                r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}= [num2cell((1:details.n_subjs)')     cell(details.n_subjs, 2)];
                
                % ######################################################
                % Expectation minimization step (one iteration of fitting all subjects)
                % ######################################################
                
                for s=1:details.n_subjs;
                    %                 disp(['   Subject ' num2str(s)])
                    wsi.i_sub=1; wsi.fitok= - 999 ;
                    
                    % Commence fit for this subject (keep going with random starts until a fit is reached, no fminunc error. This should be quick)
                    while wsi.fitok<1;
                        
                        % (a) Set parameters for fminunc, and constrain (param --> fit space)
                        %           For now, initialize subject fits with 0s. If problematic, specify with respect to sensible parameter ranges.
                        %           After 1st go, initialize fit with estimates from last known fit-params
                        if wsi.i_sub==1 && i_group==1
                            wsi.startpar=randn(1,details.models{m,2});
                        elseif wsi.i_sub==1 && i_group>1
                            wsi.startpar=r_iterations{m,1+2*(i_evmod-1)+1}{i_group-1,1}(s, 5: 4+details.models{m,2});
                        else
                            %                         wsi.parok=0;
                            %                         while wsi.parok~=0
                            if exist('wsi.fitparvals', 'var')==1
                                wsi.startpar=wsi.fitparvals+.1*randn(1,details.models{m,2});  % If this is not the 1st attempt (take the last fit and add noise)
                            else wsi.startpar=randn(1,details.models{m,2});
                            end
                            %                             wsi.checkparvals=f_transpar(details.models{m,3}, wsi.startpar, 'to');  % Are parameter values ok?
                            %                             if sum(isinf(wsi.checkparvals))==0;  wsi.parok=1; end
                            %                         end
                        end
                        % if strcmp(details.models{m,3}{p}, 'm')~=1; wsi.startpar.startpar(p)=round(wsi.startpar.startpar(p)*10)/10; end % Round
                        
                        % (b) Run nll+value function through fminunc
                        try
                            if strcmp(details.models{m,1}(2), 'p')==1
                                [wsi.fitparvals, wsi.L, wsi.fitok , wsi.output, wsi.grad , wsi.hess]=fminunc(@(x)f_nllsoftmax_lapse(x, {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 1-(i_group==1)  ), wsi.startpar, w.options);
                                %                             [nll,pch]=f_nllsoftmax_lapse(wsi.startpar, {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 1-(i_group==1)  );
                                %                             f_transpar({'b';'p';'j';'m';'f';'e';'w'},wsi.startpar, 'to')
                                %                             if rank(wsi.hess)  <  size(wsi.hess,1); keyboard; end   % singular hessian?
                            else
                                [wsi.fitparvals, wsi.L, wsi.fitok , wsi.output, wsi.grad , wsi.hess]=fminunc(@(x)f_nllsoftmax(x, {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 1-(i_group==1) ), wsi.startpar, w.options);
                                % f_nllsoftmax(wsi.fitparvals, {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 1-(i_group==1) )
                            end
                            
                            % (c) Check: has the subject-fit has converged (check fminunc exitflag)?
                            if wsi.fitok<=0 ;
                                wsi.i_sub=wsi.i_sub+1;
                                if wsi.i_sub>3; wsi.fitok=1;  errorlog{size(errorlog,1)+1,1}=['Couldn''t fit ' details.subjects{s} ' at all for  '  details.models{m,1} ' (>3 attempts)'] ; disp(errorlog{size(errorlog,1),1}); end
                            end
                        catch MExc
                            errorlog{size(errorlog,1)+1,1}=['fminunc error attepmting to fit ' details.subjects{s} ' for  '  details.models{m,1}] ; disp(errorlog{size(errorlog,1),1});
                            errorlog{end,2}= [details.models{m} ' ' details.subjects{s}];
                            errorlog{end,3}= {MExc Z wsi.startpar  i_group wsi.i_sub };
                            
                            % Have the last 1000 errors been the same? If so, just prevent memory errors
                            if size(errorlog,1)>1000 &  sum(cellfun(@(x)strcmp(x, errorlog{end,2}),  errorlog(size(errorlog,1)-999:end,2)))==1000;
                                errorlog{size(errorlog,1)+1,1}='Crashing out because the last 1000 errors have been the same';
                                try f_sendemail('kurzlich', ['[' ThisPC '] Hierarfit script crashed - last 1000 errors have been the same' ], ' ',1); end
                                %
                                error(errorlog{end,1});
                            end
                            wsi.i_sub=wsi.i_sub+1;
                        end
                    end
                    
                    % Record results of fit for this subject
                    if strcmp(details.models{m,1}(2), 'p')==1; r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 2)=f_nllsoftmax_lapse(wsi.fitparvals, {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 0);  % true nLL (no population adjustment)
                    else r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 2)=f_nllsoftmax(wsi.fitparvals, {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 0);
                    end
                    r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 3)=wsi.L; 		% posterior nLL (with population prior)
                    r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 4)=wsi.i_sub;
                    %                 r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 5: 4+details.models{m,2})=f_transpar(details.models{m,3}, wsi.fitparvals, 'to') ;
                    r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 5: 4+details.models{m,2})=wsi.fitparvals;
                    r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}{s,2}=wsi.hess;
                    r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}(s,3:2+details.models{m,2})=num2cell(diag(inv(full(wsi.hess)))');
                    wsi=[];
                end
                
                % Assemble group-level stats on parameter distributions (for this group-level iter)
                wg.pars=r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(:, 5:end);
                wg.variance=  cell2mat(r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}(:,3:2+details.models{m,2}));
                wg.hessians=r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}(:,2);
                Z.mu=mean(wg.pars);
                Z.nu=sqrt(  sum(wg.pars.^2 + wg.variance)  /details.n_subjs - mean(wg.pars) .^2   );
                Z.nui=inv(diag(Z.nu));
                r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}{details.n_subjs+1,2}=Z;    % Z parameters in fit-space
                
                % In calculating group stats, exclude subjects with bad variances (inf, negative)
                if sum(isnan(Z.nu))~=0 || sum(isinf(Z.nu))~=0 || isreal(Z.nu)==0
                    disp(['   Altering Z    ('    num2str(i_group) 'th ev step)']);  disp('     original Z:'); disp(Z)
                    wsi.subsok=ones(size(wg.pars));
                    wsi.subsok=wsi.subsok-isnan(wg.variance)-(wg.variance<0)-isinf(wg.variance); % disqualifying criteria
                    wsi.subsok=wsi.subsok==ones(size(wg.pars));
                    for k=1:details.models{m,2}
                        wk.v=wg.variance(find(wsi.subsok(:,k)), k); wk.pars=wg.pars(find(wsi.subsok(:,k)), k);
                        Z.mu(k)=mean( wk.pars );
                        Z.nu(k)=sqrt(  sum(wk.pars .^2 + wk.v)  /length(wk.pars) - mean(wk.pars) .^2   );
                        wk=[];
                    end
                    Z.nui=inv(diag(Z.nu));
                    disp('     NewZ:'); disp(Z)
                end
                
                % Assemble group-level fit (Ls, mus and nus)
                %   (1)= sum of nLL actuals, (2)=sum of nLL posterior, 3 onwards:  mu's and nu's for each param (all mus, then all nus)
                we.iters(i_group, :)= [sum(r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(:, 2)) sum(r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(:, 3)) Z.mu Z.nu];
                r_iterd{m,2+i_evmod}(i_group, 1:length(we.iters(i_group, :)))=we.iters(i_group, :);
                
                if  i_group>1 && sum(   abs(  we.iters(i_group, 2:end) - we.iters(i_group-1, 2:end) )   )   <1e-3;
                    disp(['   [' num2str(m) ' of '  num2str(length(details.whichmodels)) '] Group converged on ' num2str(i_group) 'th try         (mu:  ' num2str(Z.mu) '       nu:  ' num2str(Z.nu)  ')']);
                    details.modelgroupfitok(m)=1;
                    break;
                else i_group=i_group+1;
                    if mod(i_group, 10)==0;   disp(['         Attempting group fit (within EV): '  num2str(i_group) 'th try         (mu:  ' num2str(Z.mu,2) '       nu:  ' num2str(Z.nu,2)  ')']);  end
                    if i_group>details.maxgroupiter;
                        details.modelgroupfitok(m)=-999; we=[];
                        errorlog{size(errorlog,1)+1} =[ details.whichmodels{m} ':  >' num2str(details.maxgroupiter) ' group-level iters couldn''t converge'];
                        try f_sendemail('kurzlich', ['[' ThisPC '] Modelling hierarfit: model cannot converge (' details.models{m,1} ')' ], ' ',1); end
                        break
                    end
                end
                
                % Reset Z? Decide and implement
                if sum(isnan(Z.nu(:)))~=0+sum(isinf(Z.nu(:)))~=0+sum(isreal(Z.nu(:))==0)~=0+sum(Z.nu(:)<0)~=0    || sum(Z.nu>1000000000000)>0
                    wg.resetZ=1; % Bad Zs
                elseif sum(cell2mat(cellfun(@(x)sum(isnan(x(:))),wg.hessians(1:details.n_subjs), 'UniformOutput',0))) ~=0  ||   sum(cell2mat(cellfun(@(x)sum(isinf(x(:))),wg.hessians(1:details.n_subjs), 'UniformOutput',0))) ~=0
                    % Individual subjects have bad hessians. THINK THIS THROUGH. Maybe Zs shouldn't be reset just for individual subjects ???
                    % Maybe this is why some models are failing to fit entirely
                    wg.resetZ=1;
                elseif  sum(cellfun(@(x)rank(x),wg.hessians(1:details.n_subjs))<details.models{m,2})==details.n_subjs  % All hessians are rank
                    wg.resetZ=1;
                else wg.resetZ=0;
                end
                if wg.resetZ==1;  % Implement reset
                    i_group=i_group-1;
                    details.resetZ{size(details.resetZ,1)+1, 1}=details.models{m,1};
                    details.resetZ{end, 2}=wg;
                    disp('Resetting Z. Hessians rank for all subjects, or variances are bad even after correction')
                    if mod(sum(strcmp(details.resetZ(:,1), details.models{m,1})), 10)==0; try w.c=clock; f_sendemail('kurzlich', ['[' ThisPC '] Hierarchical fit reset for model ' details.models{m,1} '  (reset '  num2str(sum(strcmp(details.resetZ(:,1), details.models{m,1}))) 'x)  [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end; end
                    
                    % If Z has been serially reset too many times, break
                    if size(details.resetZ,1)>100 &&  sum(cellfun(@(x)strcmp(x, details.resetZ{end,1}),  details.resetZ(size(details.resetZ,1)-49:end,1)))==0;
                        details.modelgroupfitok(m)=-999;
                        errorlog{size(errorlog,1)+1,1} = [ details.whichmodels{m} ':  Group-level fit aborted because Z had to be reset too many times. The model is likely problematic!']; disp(errorlog{end,1})
                        errorlog{end, 2} =[ details.whichmodels{m} ' Z reset too many times'];
                        errorlog{end, 3} =Z;
                    end
                else wg=[];
                end
                
                
                % Give up on this model fit?
                if details.modelgroupfitok(m)==-999;  break; end
            end
            
            
            
            % ######################################################
            % At the end of ONE successful group-level EV fit, calculate BIC for the model
            % ######################################################
            disp('   Calculating BIC via samples'); rand('seed',sum(100*clock));
            if details.modelgroupfitok(m)==1
                for s=1:details.n_subjs
                    ws.sampleparvals = diag(sqrt(Z.nu)) * randn(details.models{m,2},  details.calcbic_n_samples)   +  Z.mu' * ones(1, details.calcbic_n_samples);    % Sampled params for this subject, given the group distribution
                    ws.sampleparvals=sortrows(ws.sampleparvals',1)';
                    ws.LLi=zeros(details.calcbic_n_samples,1);
                    
                    for k=1:details.calcbic_n_samples
                        if strcmp(details.models{m,1}(2), 'p')==1  % Calculate nLL for sampled subject-parameters
                            ws.LLi(k)=f_nllsoftmax_lapse(ws.sampleparvals(:, k)', {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 0);  % true nLL (no population adjustment)
                        else ws.LLi(k)=f_nllsoftmax(ws.sampleparvals(:, k)', {details.models{m,1} subjdata{s,1+details.tasktype} details.fixedpar col},Z, 0);
                        end
                    end
                    ws.LLi=ws.LLi(isinf(ws.LLi)==0 & isnan(ws.LLi)==0);   % Remove bad nLLs (infinites & nans)
                    if length(ws.LLi)< details.calcbic_n_samples*0.01; disp(['      >10% bad samples removed in calculating BIC via sampling (' details.subjects{s} ')']); end
                    
                    we.iL(s) = log(sum(exp(-ws.LLi))/details.calcbic_n_samples);  % Mean nLL for this subject, for this parameter distribution
                    ws=[];
                end
                
                % End of 1 EV-iteration for entire group x model!!! ----------------------------------
                r_iterations{m,1+2*(i_evmod-1)+2}= -2*sum(we.iL)   + details.models{m,2}*log(sum(details.subj_ntrials));
            end
            we=[];
            
        end
        %   i_group/wg= iteration within EV process (ev step, automated),  i_evmod/we=requested iteration of entire model fit (outside group EV protocol)
        if details.modelgroupfitok(m)==1
            % ######################################################
            % Pick best iteration from requested iterations (lowest BIC)
            % ######################################################
            wm.best=find(cell2mat(r_iterations(m,3:2:size(r_iterations,2)))==min(cell2mat(r_iterations(m,3:2:size(r_iterations,2)))));
            r_res{m,2}(1:details.n_subjs, 1)=1:details.n_subjs;
            r_res{m,2}(:, 2)=r_iterations{m,1+2*(wm.best-1)+1}{end,1}(:, 2); % recorded once: correct L
            r_res{m,2}(:, 4:3+details.models{m,2})=r_iterations{m,1+2*(wm.best-1)+1}{end,1}(:, 5:end);
            r_res{m,3}=r_iterations{m,1+2*(wm.best-1)+2};
            r_res{m,4}=r_iterations{m,1+2*(wm.best-1)+1}{end,2}(:, 2);
            r_res(m,6:5+details.models{m,2})=num2cell(mean(r_res{m,2}(:, 4:3+details.models{m,2})));
            wm=[];
        end
        
        % Save partial fits
        w.c=clock; save(['1 Inputs' filesep '1b PartialFits' filesep 'Partial hierarfit ' details.task ' (' date ' - ' num2str(w.c(4)) ' hrs)']);
    end
    try w.c=clock; f_sendemail('kurzlich', ['[' ThisPC '] Slack hierarfit almost done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end
    
end

if details.do.doreturn
    
    % Return all parameters (iterations) into true parameter space
    disp('Returning parameters from fit-space to true paramenter space ##############')
    for m=startfrom:details.n_models
        if details.modelgroupfitok(m)==1
            % r_iteration
            for i_evmod=1:details.fitgroup_n_eviter
                for i_group=1:size(r_iterations{m,1+2*(i_evmod-1)+1},1)
                    
                    % Correct parameters
                    for s=1:details.n_subjs
                        r_iterations{m, 1+2*(i_evmod-1)+1}{i_group,1}(s, 5:end)=f_transpar(details.models{m,3}, r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(s, 5:end), 'to');
                    end
                    
                    % Characterize distribution
                    r_iterations{m,1+2*(i_evmod-1)+1}{i_group,2}{details.n_subjs+1,2}.parspace_mu=mean(r_iterations{m,1+2*(i_evmod-1)+1}{i_group,1}(:, 5:end));
                end
            end
            
            % r_res
            for s=1:details.n_subjs
                r_res{m,2}(s,4:end)=f_transpar(details.models{m,3}, r_res{m,2}(s,4:end), 'to');
            end
        end
    end
    for m=1:size(r_res,1) % Mean parameter values in real space
        if details.modelgroupfitok(m)==1
            for p=1:size(r_res{m,2},2)-3
                r_res{m, 5+p}=mean(r_res{m,2}(:, p+3));
            end
        else r_res{m, 3}=nan;
        end
    end
    r_res(find(cellfun(@(x)isempty(x), r_res(:,3))), 3)=num2cell(repmat(-999, sum(cellfun(@(x)isempty(x), r_res(:,3))), 1));
    r_res=sortrows(r_res,3);
    
    % Collate data for plotting convergence
    for m=startfrom:details.n_models
        if details.modelgroupfitok(m)==1
            for p=1:details.models{m,2}
                r_iterd{m,2}{p}=nan(details.n_subjs, size(r_iterd{m,3},1));
                for i=1:size(r_iterd{m,3},1)
                    r_iterd{m,2}{p}(:,i)=r_iterations{m,2}{i,1}(:, p+4);
                end
            end
        end
    end
    
end

%% END 

disp('===================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp('INSTRUCTION: Results variables (prefix r):r_iterations, r_res, r_models')
disp('                     [r_res: Col 1= Model, Col 2=Subject fits, 3=Model BIC, Col 4=N valid subjects, Col 6 onwards=Mean parameter values]'); disp(' ')
disp('Check command-window log for poor fits; See ''further instructions'' at end of script for Bayes factor & model weights'); disp(' ')
disp('====================================')

% Save
resfilename=['res_hierarfitmodels_' details.task ' (' date ')'];
resfilewhere=[where.where filesep '1 Inputs' filesep details.folder4saving  filesep];  cd(resfilewhere), filelist=cellstr(ls);
if sum(strcmp(filelist, [resfilename '.mat']))>0;  k=2; kk=0;
    while kk==0
        if sum(strcmp(filelist, [resfilename  num2str(k) '.mat']))==0
            resfilename=[resfilename num2str(k)]; kk=1;
        else k=k+1;
        end
    end
end
save([resfilewhere resfilename], 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog'); diary off
try % Transfer file to DeletedDaily 
    movetowhere='\\Asia\DeletedDaily\EL ModRuns';
    if isdir(movetowhere)==0; mkdir(movetowhere);end
    w.c=clock; cd(resfilewhere); 
    transferfilename=[ThisPC '_' num2str(w.c(3)) '_' num2str(w.c(2)) 'at' num2str(w.c(4)) num2str(w.c(5)) 'hrs ' resfilename];
    copyfile([resfilename '.mat'],  [movetowhere filesep transferfilename '.mat']);        
    disp(['Saved locally:  '  resfilename ]); disp(['Saved in DeletedDaily:  '  transferfilename])
catch
    disp(['Saved locally:  '  resfilename ' - NOT transferred to DeletedDaily']);
end


%% Further instructions (Copy to command window to calculate

openvar r_res; diary off; char(errorlog{:,1})
disp('Any models failed to converge? ######');  disp([details.whichmodels num2cell(details.modelgroupfitok)])


% Plot convergence for best model
modname=r_res{1,1};
% modname='b16_feow';
m=find(strcmp(r_iterd(:,1),  modname)); 
if isempty(m)==1; error('Requested invalid model name for plotting!'); end; figure('color','w');
for p=1:  details.models{m,2}
    subplot(1, details.models{m,2}, p)
    plot(r_iterd{m,2}{p}')
    title([r_iterd{m,1} ' - '   details.models{m,3}{p}])
end

% Calculate Bayes factor ---------------------- 
%       Which models to compare? (Row number in 'r_res')
m1=1;
m2=2;
B=(r_res{m1,3}-r_res{m2,3})*-0.5;
L= sum(r_res{1,2}(:,2));   R = sum(details.subj_ntrials) .* -log(1/3);

disp(['B=' num2str(B)  '  (m1=' r_res{m1,1} ', m2='  r_res{m2,1} ')'])
disp(['Pseudo r2 for winning model= :' num2str(   1-L/R )]);
% Interpreting B (conventions): 3-10=moderate evidence, >10=strong evidence (in favour of m1)


