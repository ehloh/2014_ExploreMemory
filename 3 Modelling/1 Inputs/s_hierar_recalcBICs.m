% Instructions: Manually load the hierarchical file
clear all; clc

root='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling';

% root='/Users/EleanorL/Dropbox/sandisk/8 explore mem/3b modelling';


cd([ root '\1 Inputs\2b Hierar fits\New\v2'])
mkdir('BIC_PartialFits')

filename='res_hierarfitmodels_cF (16-Apr-2015) v2 remaining';
datafile='All data v2 (01-Apr-2015)';

nsamples=10000;

%%
load([filename '.mat'])
s=load([root filesep '1 Inputs' filesep datafile '.mat']); subjdata=s.subjdata; col=s.details.col;
valwhere=[root filesep '2 Val fxn']; modfol=root; 
path(pathdef), addpath(valwhere);  addpath(genpath(valwhere)); addpath(modfol)
 

% 
details.models=sortrows(details.models,1);  r_res=sortrows(r_res,1); r_iterations=sortrows(r_iterations,1); r_iterd=sortrows(r_iterd,1); details.whichmodels=sortrows(details.whichmodels,1);
[ r_res , r_iterations, r_iterd, missingres] = f_hierarmatch(r_res , r_iterations, r_iterd);   % Check that material for resampling is there
if isempty(missingres)==0; openvar missingres, error('Some mismatch between r_res , r_iterations, r_iterd! Check files first (see missingres)'); end
%
details.calcbic_n_samples=nsamples;

for m=1:size(r_res,1)
    
    % Plug in details from model fit process
    numus= r_iterd{  strcmp(r_iterd(:,1), r_res{m,1}) ,3}(end, 3:end);  % (1)= sum of nLL actuals, (2)=sum of nLL posterior, 3 onwards:  mu's and nu's for each param (all mus, then all nus)
    Z.mu = numus(1:length(numus)/2);
    Z.nu = numus(length(numus)/2+1:end);
    i_evmod=1;
    
    
    
    %% Section from within-EV code
    
    disp([ num2str(m) ' out of ' num2str(size(r_res,1)) ':  Calculating BIC via samples for ' r_res{m,1}]); rand('seed',sum(100*clock));
    for s=1:  details.n_subjs
        disp(s)
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
    bic= -2*sum(we.iL)   + details.models{m,2}*log(sum(details.subj_ntrials));
    r_iterations{strcmp(r_iterations(:,1), r_res{m,1}),1+2*(i_evmod-1)+2}=bic;
    
    %% Plug re-calculated BICs into the rest of the code
    r_res{m,3}=bic; 
    
    % SAVE partialfits
    w.c=clock; save(['BIC_PartialFits' filesep filename ' Partial (' date ' - ' num2str(w.c(4)) ' hrs)'])
    
    
    we=[];
end


% SAVE
try w.c=clock; f_sendemail('kurzlich', ['Hierarchical recalc BIC done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end
resfilename=f_newname([filename ' BICredone'], pwd);
disp(['Name of new file:   ' resfilename]), disp(['Directory  : ' pwd])
input('SAVE?   ');
save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog'); diary off
disp('Done - don''t forget to move over to correct folder!')


