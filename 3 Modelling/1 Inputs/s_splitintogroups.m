% Combine runs (fitmodels, not hierarfits) that were run separately in execution
clear all, close all hidden, clc

% REQUEST FILE
request.resfilename='2a Fixed fits\res_fitmodels_cF (05-Apr-2015) all bpjm16_L6472';
%
request.vok{1}={'p01_AF' ; 'p01_AS' ; 'p02_TY' ; 'p03_JV' ; 'p06_MK' ; 'p08_AK' ; 'p10_SG' ; 'p11_OG' ; 'p30_AR' ; 'p31_RD' ; 'p33_AA' ; 'p34_BE'};
request.vok{2}={'p14_CP';'p15_SZ';'p23_SG';'p24_VP';'p26_PB';'p27_PV';'p29_MY';'p35_EO';'p37_ND';'p38_SL';'p39_JH';'p40_CC'; 'p45_SO';};

%%

% 
r=load([request.resfilename '.mat']); newr=cell(2,1);
for v=1:2 % Get subject no.s
    request.subno{v}=cellfun(@(x)find(strcmp(r.details.subjects, x)), request.vok{v});
end
for v=1:2
    for m=1:size(r.r_res,1)
        % 'r_iterations' - model fits for each subject x iteration
        %       Col 1=Model name
        %       Col 2=Iteration results (nLL, parameters)
        %             Col 1: 
        %             Col 2: nLL
        %             Col 3: Iterations exceeded fminunc?
        %             Col 4 onwards: model parameters (1st parameter is beta/inverse temperature)
        %       Col 3=Iteration hessians
        %             Col 2: Hessian for this iteration
        %       Col 4=nLL histogram data
        %       Col 5= Subject hessians for best fit (col 1), ok? (non-singular; col 2)
        newr{v}.r_iterations{m,1}=r.r_iterations{m,1};
        newr{v}.r_iterations{m,2}=[]; newr{v}.r_iterations{m,3}=[]; newr{v}.r_iterations{m,4}=[];
        newr{v}.r_iterations{m,5}=r.r_iterations{m,5}(request.subno{v},:);
        for s=1:length(request.subno{v})
            newr{v}.r_iterations{m,2}=[newr{v}.r_iterations{m,2}; r.r_iterations{m,2}(r.r_iterations{m,2}(:,1)==request.subno{v}(s),:)];
            newr{v}.r_iterations{m,3}=[newr{v}.r_iterations{m,3}; r.r_iterations{m,3}(cell2mat(r.r_iterations{m,3}(:,1))==request.subno{v}(s),:)];
            % Col 4 left out 
        end
        
        % 'r_res'
        %       Col 1: Model name
        %       Col 2: Subject fit parameters
        %             Col i: BIC
        %             Col ii: nLL
        %             Col iii: fminunc exceeded default iterations?
        %             Col iv onwards: parameters (beta first)
        %       Col 3: Model BIC (summed across subjects)
        %       Col 4: Hessians
        %       Col 5: 
        %       Col 6 onwards: mean parameter values
        newr{v}.r_res{m,1}=r.r_res{m,1};
        newr{v}.r_res{m,2}=r.r_res{m,2}(request.subno{v},:);
        newr{v}.r_res{m,4}=r.r_res{m,4}(request.subno{v},:);
        
        % Re-calculate group stats per model
        newr{v}.r_res{m,3}=sum(newr{v}.r_res{m,2}(:,1));
        for p=1: size(newr{v}.r_res{m,2},2)-3
            newr{v}.r_res{m,5+p}=mean(newr{v}.r_res{m,2}(:, 3+p));
        end
    end
    newr{v}.r_res=sortrows(newr{v}.r_res,3);
    
    % Compile other details
    newr{v}.details=r.details;
    newr{v}.details.subjects=r.details.subjects(request.subno{v},:);
    newr{v}.details.n_subjs=length(newr{v}.details.subjects);
    newr{v}.details.dataset=['Split into groups from: ' newr{v}.details.dataset];
    
    % Save
    input(['SAVE version ' num2str(v) '?']);
    details=newr{v}.details; r_res=newr{v}.r_res; r_iterations=newr{v}.r_iterations; errorlog=r.errorlog; rc=r.rc;
    save([request.resfilename ' v' num2str(v)], 'details','r_res','r_iterations','errorlog')
end
