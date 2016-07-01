% Combine runs (fitmodels, not hierarfits) that were run separately in execution
clear all, close all hidden, clc

% [ INSTRUCTIONS ] ####################
%   1: Turn on dofit/dohierarfit
%   2: Specify which fits to combine (in the fit/hierarfit module)
%
% ---------------------------------------------------------------

dofit=0; 
dohierarfit=1;

%% Combine runs from fminunc fit

% Request #####
log.runs2combine={
        %     'res_fitmodels_cF (04-Apr-2015) bjmi bpjmi'
        %     'res_fitmodels_cF (04-Apr-2015) bpji'
        %     'res_fitmodels_cF (05-Apr-2015) bji'
    };

if dofit
    for o1=1:1
    
    for o=1:1  % General setup
        
        
        % Append full paths
        where.where='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling';
%         where.where='/Users/EleanorL/Dropbox/SANDISK/8 Explore Mem/3b Modelling';
        path(pathdef); addpath(where.where)
%         log.runs2combine=cellfun(@(x)[where.where filesep x], log.runs2combine, 'UniformOutput',0);
        log.runs2combine=cellfun(@(x)[where.where filesep '1 Inputs' filesep '2a Fixed fits' filesep x], log.runs2combine, 'UniformOutput',0);
        
        
        % Check compatible inputs
        if length(log.runs2combine)~= length(cell2mat(strfind(log.runs2combine, 'fitmodels_'))), error('Wrong type of fits requested'); end
        if isempty(strfind(log.runs2combine{1}, 'fitmodels_cF'))==0
            if length(log.runs2combine)~= length(cell2mat(strfind(log.runs2combine, 'fitmodels_cF'))), error('NOT all cF models!!!'); end
        elseif length(log.runs2combine)~= length(cell2mat(strfind(log.runs2combine, 'fitmodels_ct'))), error('NOT all ct models!!!');
        end
    end
    
    % Combine, including all relevant information (If particular models are repeated, they are just entered twice)
    d_newfit.details.fitsfrom=log.runs2combine; d_newfit.details.ReadMe={};
    d_newfit.r_res={}; d_newfit.r_iterations=[]; d_newfit.errorlog=[]; d_newfit.rc=[];
    d_allfits=cell(length(log.runs2combine),1);
    for f=1:length(log.runs2combine)  % Combine runs
        d_allfits{f}=load(log.runs2combine{f});
        d_newfit.details.fitsfrom{f,2}=d_allfits{f}.details;
        
        % Append r_res
        wf.ncols=size(d_allfits{f}.r_res,2);
        wf.startrows=size(d_newfit.r_res,1)+1;
        wf.endrows=size(d_newfit.r_res,1)+size(d_allfits{f}.r_res,1);
        d_newfit.r_res(wf.startrows:wf.endrows, 1:wf.ncols)  = d_allfits{f}.r_res;
        try d_newfit.r_iterations=[d_newfit.r_iterations;  d_allfits{f}.r_iterations];  % Append r_iterations
        catch d_newfit.details.ReadMe{size(d_newfit.details.ReadMe,1)+1}=['r_iterations missing for fit #' num2str(f) '(see details.fitsfrom)'];
        end
        try d_newfit.errorlog=[d_newfit.errorlog; d_allfits{f}.errorlog];
        catch d_newfit.details.ReadMe{size(d_newfit.details.ReadMe,1)+1}=['errorlog missing for fit #' num2str(f) '(see details.fitsfrom)'];
        end
        
        
    end
    for o=1:1 % Record details
        d_newfit.details.tasktype=d_allfits{1}.details.tasktype;% Details that assumed to be the same
        d_newfit.details.task=d_allfits{1}.details.task;
        d_newfit.details.modelfams=d_allfits{1}.details.modelfams;
        d_newfit.details.model_defaults=d_allfits{1}.details.model_defaults;
        d_newfit.details.par_transformations=d_allfits{1}.details.par_transformations;
        d_newfit.details.fixedpar=d_allfits{1}.details.fixedpar;
        d_newfit.details.col=d_allfits{1}.details.col;
        try d_newfit.rc=d_allfits{1}.details.rc;  catch d_newfit.rc='Unspecified'; end
        
        % Details that need to be verified
        if length(unique(cellfun(@(x)x.details.n_iterations, d_allfits)))==1; d_newfit.details.n_iterations=d_allfits{1}.details.n_iterations; else  d_newfit.details.n_iterations='See individual fits (details.fitsfrom)'; end
        if length(unique(cellfun(@(x)x.details.n_subjs, d_allfits)))==1;  d_newfit.details.n_subjs=d_allfits{1}.details.n_subjs; else  d_newfit.details.n_subjs='See individual fits (details.fitsfrom)'; end
        if length(unique(cell2mat(cellfun(@(x)length(x.details.subjects), d_allfits, 'UniformOutput',0))))==1;   d_newfit.details.subjects=d_allfits{1}.details.subjects; else  d_newfit.details.subjects='See individual fits (details.fitsfrom)'; end
        if sum(strcmp(cellfun(@(x)(x.details.dataset), d_allfits, 'UniformOutput',0), d_allfits{1}.details.dataset))==length(log.runs2combine);  d_newfit.details.dataset=d_allfits{1}.details.dataset;  else  input('Different datasets! you sure you want to combine?'); d_newfit.details.dataset='See individual fits (details.fitsfrom)'; end
        if sum(strcmp(cellfun(@(x)(x.details.dataset_date), d_allfits, 'UniformOutput',0), d_allfits{1}.details.dataset_date))==length(log.runs2combine);  d_newfit.details.dataset_date=d_allfits{1}.details.dataset_date;  else  d_newfit.details.dataset_date='See individual fits (details.fitsfrom)'; end
%         if sum(strcmp(cellfun(@(x)(x.details.Vanfxn_type), d_allfits, 'UniformOutput',0), d_allfits{1}.details.Vanfxn_type))==length(log.runs2combine); d_newfit.details.Valfxn_type=d_allfits{1}.details.Vanfxn_type;  else  d_newfit.details.Valfxn_type='See individual fits (details.fitsfrom)'; end
% if sum(strcmp(cellfun(@(x)(x.details.Valfxn_type), d_allfits, 'UniformOutput',0), d_allfits{1}.details.Valfxn_type))==length(log.runs2combine); 
% d_newfit.details.Valfxn_type=d_allfits{1}.details.Valfxn_type;  
% else  d_newfit.details.Valfxn_type='See individual fits (details.fitsfrom)'; 
% end
        if sum(strcmp(cellfun(@(x)(x.details.folder4saving), d_allfits, 'UniformOutput',0), d_allfits{1}.details.folder4saving))==length(log.runs2combine);  d_newfit.details.folder4saving=d_allfits{1}.details.folder4saving;  else  d_newfit.details.folder4saving='See individual fits (details.fitsfrom)'; end
        
        d_newfit.details.Valfxn_type=d_allfits{1}.details.Valfxn_type;  
        
        % Details to recompile from new consolidate4d fits:
        d_newfit.details.n_models= size(d_newfit.r_res,1);
        d_newfit.details.whichmodels=sortrows(d_newfit.r_res(:,1));
        d_newfit.details.models(:,1)=d_newfit.details.whichmodels;
        [wf.model_defaults  wf.par_transformations wf.models] = f_modelsettings(d_newfit.details.whichmodels, 20);  % Arbitrary n iterations
        d_newfit.details.models=wf.models;
    end
    
    % Outputs [Save in folder of 1st fit file]
    where.save=log.runs2combine{1}(1:cellfun(@(x)x(end), {strfind(log.runs2combine{1}, filesep)}));   
    name2save=f_newname(['res_fitmodels_' d_newfit.details.task ' (' date ')' num2str(k)], where.save);
    %
    clear('details','errorlog','r_res','r_iterations', 'rc'); details=d_newfit.details;  r_res=sortrows(d_newfit.r_res,3);   r_iterations=d_newfit.r_iterations; errorlog=d_newfit.errorlog; rc=d_newfit.rc;
    disp('Save at:'); disp(where.save ),  disp(['File name:  ' name2save]); input('Continue to save?');
    save([where.save name2save],  'details','errorlog','r_res','r_iterations', 'rc');
    end
end

%% Combine runs for hierarfit


% Request #####
log.runs2combine={
    
'2b Hierar fits\New\res_hierarfitmodels_cF (16-Apr-2015) v1 remaining BICredone'

% '2b Hierar fits\res_hierarfitmodels_cF (13-Apr-2015) v2 bpji7e_L3523 BICredone';


'2b Hierar fits\res_hierarfitmodels_cF (13-Apr-2015) v1 bpji16_L3090 BICredone'

};

if dohierarfit
    
    for o=1:1  % General setup
        % Append full paths
        where.where='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling'; where.origmod='D:\Dropbox\SANDISK\4 Explore experiment\3 Analysis\4 Fit computational models';
%         where.where='/Users/EleanorL/Dropbox/SANDISK/8 Explore Mem/3b Modelling'; where.origmod='/Users/EleanorL/Dropbox/SANDISK/4 Explore experiment/3 Analysis/4 Fit computational models';
        path(pathdef);  addpath(where.where)
        log.runs2combine=cellfun(@(x)[where.where filesep '1 Inputs' filesep x '.mat'], log.runs2combine, 'UniformOutput',0);
        
        % Check compatible inputs
        if length(log.runs2combine)~= length(cell2mat(strfind(log.runs2combine, 'hierarfitmodels_'))), error('Wrong type of fits requested!'); end
        if isempty(strfind(log.runs2combine{1}, 'hierarfitmodels_cF'))==0
            if length(log.runs2combine)~= length(cell2mat(strfind(log.runs2combine, 'hierarfitmodels_cF'))), error('NOT all cF models!!!'); end
        elseif length(log.runs2combine)~= length(cell2mat(strfind(log.runs2combine, 'hierarfitmodels_ct'))), error('NOT all ct models!!!');
        end
    end
    
    % Combine, including all relevant information (If particular models are repeated, they are just entered twice)
    d_newfit.details.fitsfrom=log.runs2combine;
    d_newfit.r_res={}; d_newfit.r_iterd=[];  d_newfit.r_iterations=[]; d_newfit.errorlog=[]; d_newfit.rc=[];
    d_allfits=cell(length(log.runs2combine),1);
    for f=1:length(log.runs2combine)  % Combine runs
        d_allfits{f}=load(log.runs2combine{f});
        d_newfit.details.fitsfrom{f,2}=d_allfits{f}.details;
        
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
        wf.ncols=size(d_allfits{f}.r_res,2);
        wf.startrows=size(d_newfit.r_res,1)+1;
        wf.endrows=size(d_newfit.r_res,1)+size(d_allfits{f}.r_res,1);
        d_newfit.r_res(wf.startrows:wf.endrows, 1:wf.ncols)  = d_allfits{f}.r_res;
        
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
        %       Col 2               Plot matrix (row=subject, col=iteration)
        %       Col 3               Details for requested iter 1..
        d_newfit.r_iterations=[d_newfit.r_iterations;  d_allfits{f}.r_iterations];
        d_newfit.r_iterd=[d_newfit.r_iterd;  d_allfits{f}.r_iterd];
        
        
        
        if f>1 & size(d_newfit.errorlog,2) < size(d_allfits{f}.errorlog,2)
            d_newfit.errorlog{end, size(d_allfits{f}.errorlog,2)}=1;
        elseif  f>1 & size(d_newfit.errorlog,2) > size(d_allfits{f}.errorlog,2)
%             if f==2, 
%                 d_allfits{f}.errorlog{1, size(d_newfit.errorlog,2)}=1; 
% 
%             else
                                d_allfits{f}.errorlog{end, size(d_newfit.errorlog,2)}=1;
%             end
        end
        d_newfit.errorlog=[d_newfit.errorlog; d_allfits{f}.errorlog];
    end
    for o=1:1 % Record details
        d_newfit.details.tasktype=d_allfits{1}.details.tasktype;% Details that assumed to be the same
        d_newfit.details.task=d_allfits{1}.details.task;
        d_newfit.details.modelfams=d_allfits{1}.details.modelfams;
        d_newfit.details.model_defaults=d_allfits{1}.details.model_defaults;
        d_newfit.details.par_transformations=d_allfits{1}.details.par_transformations;
        d_newfit.details.fixedpar=d_allfits{1}.details.fixedpar;
        d_newfit.details.col=d_allfits{1}.details.col;
        d_newfit.details.resetZ=d_allfits{1}.details.resetZ;
        d_newfit.details.subj_ntrials=d_allfits{1}.details.subj_ntrials;
        try d_newfit.rc=d_allfits{1}.details.rc;  catch d_newfit.rc='Unspecified'; end
        
        % Details that need to be verified
        if length(unique(cellfun(@(x)x.details.n_subjs, d_allfits)))==1;  d_newfit.details.n_subjs=d_allfits{1}.details.n_subjs; else  d_newfit.details.n_subjs='See individual fits (details.fitsfrom)'; end
        if length(unique(cell2mat(cellfun(@(x)length(x.details.subjects), d_allfits, 'UniformOutput',0))))==1;   d_newfit.details.subjects=d_allfits{1}.details.subjects; else  d_newfit.details.subjects='See individual fits (details.fitsfrom)'; end
        if sum(strcmp(cellfun(@(x)(x.details.dataset), d_allfits, 'UniformOutput',0), d_allfits{1}.details.dataset))==length(log.runs2combine);  d_newfit.details.dataset=d_allfits{1}.details.dataset;  else  input('Different datasets! you sure you want to combine?'); d_newfit.details.dataset='See individual fits (details.fitsfrom)'; end
        if sum(strcmp(cellfun(@(x)(x.details.dataset_date), d_allfits, 'UniformOutput',0), d_allfits{1}.details.dataset_date))==length(log.runs2combine);  d_newfit.details.dataset_date=d_allfits{1}.details.dataset_date;  else  d_newfit.details.dataset_date='See individual fits (details.fitsfrom)'; end
        %         if sum(strcmp(cellfun(@(x)(x.details.Valfxn_type), d_allfits, 'UniformOutput',0), d_allfits{1}.details.Valfxn_type))==length(log.runs2combine); d_newfit.details.Valfxn_type=d_allfits{1}.details.Vanfxn_type;  else  d_newfit.details.Valfxn_type='See individual fits (details.fitsfrom)'; end
        %         if length(unique(cellfun(@(x)x.details.n_iterations, d_allfits)))==1; d_newfit.details.n_iterations=d_allfits{1}.details.n_iterations; else  d_newfit.details.n_iterations='See individual fits (details.fitsfrom)'; end
        if sum(strcmp(cellfun(@(x)(x.details.folder4saving), d_allfits, 'UniformOutput',0), d_allfits{1}.details.folder4saving))==length(log.runs2combine);  d_newfit.details.folder4saving=d_allfits{1}.details.folder4saving;  else  d_newfit.details.folder4saving='See individual fits (details.fitsfrom)'; end
        if length(unique(cellfun(@(x)x.details.fitgroup_n_eviter, d_allfits))) ==1;   d_newfit.details.fitgroup_n_eviter=d_allfits{1}.details.fitgroup_n_eviter;  else  d_newfit.details.fitgroup_n_eviter='See individual fits (details.fitsfrom)'; end
        if length(unique(cellfun(@(x)x.details.calcbic_n_samples, d_allfits))) ==1;   d_newfit.details.calcbic_n_samples=d_allfits{1}.details.calcbic_n_samples;  else  d_newfit.details.calcbic_n_samples='See individual fits (details.fitsfrom)'; end
        if length(unique(cellfun(@(x)x.details.fitsub_n_iter, d_allfits))) ==1;   d_newfit.details.fitsub_n_iter=d_allfits{1}.details.fitsub_n_iter;  else  d_newfit.details.fitsub_n_iter='See individual fits (details.fitsfrom)'; end
        
        % Details to recompile from new consolidate4d fits:
        d_newfit.details.n_models= size(d_newfit.r_res,1);
        d_newfit.details.whichmodels=sortrows(d_newfit.r_res(:,1));
        %         d_newfit.details.models(:,1)=d_newfit.details.whichmodels;
        [wf.model_defaults  wf.par_transformations wf.models] = f_modelsettings(d_newfit.details.whichmodels, 20);  % Arbitrary n iterations
        d_newfit.details.models=wf.models;
    end
    
    % Outputs [Save in folder of 1st fit file]
    where.save=log.runs2combine{1}(1:cellfun(@(x)x(end), {strfind(log.runs2combine{1}, filesep)}));  
    name2save=f_newname(['res_hierarfitmodels_' d_newfit.details.task ' (' date ')'], where.save);
    clear('details','errorlog','r_res','r_iterations', 'r_iterd', 'rc');
    r_iterd=d_newfit.r_iterd; details=d_newfit.details;   r_res=d_newfit.r_res;
%     r_res=sortrows(r_res,3);   
    r_iterations=d_newfit.r_iterations; errorlog=d_newfit.errorlog; rc=d_newfit.rc;
    disp('Save at:'); disp(where.save ),  disp(['File name:  ' name2save]); input('Continue to save?');
    save([where.save name2save],  'details','errorlog','r_res','r_iterations', 'r_iterd','rc');
end

% What are the datasets?
   cellfun(@(x)(x.details.dataset), d_allfits, 'UniformOutput',0)
   