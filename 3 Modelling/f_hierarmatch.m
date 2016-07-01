function [ nrr, nri nrd, missingres] = f_hierarmatch(rr, ri, rd)
% [ r_res , r_iterations, r_iterd, missingres] = f_hierarmatch(r_res , r_iterations, r_iterd)
% Make sure that the model results listed in r_res match up with r_iterd
% and r_iterations. Hierarfits only. Output results include only models
% that are converged and that have their corresponding workspaces (r_iterations, r_iterd).
% missingres lists models that need patching (find workspaces), or that did
% not converge at all
%
% [ASSUMED ORGANIZATION] ############
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
%
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
%
% r_iterd - details of each ev step fit (within each requested ev procedure)
%       Col 1               Model name
%       Col 2               Data for plotting convergence (row=subject, col=iteration)
%       Col 3               Details for requested iter 1..

missingres={};  nrr={}; nri={}; nrd={}; modnum=1; 
for m=1:size(rr,1)
    wm.bic=rr{m,3};
    wm.modname=rr{m,1};
    if isempty(wm.bic)==0
        wm.ri_num=find( strcmp(ri(:,1), wm.modname).* cell2mat(cellfun(@(x) isempty(x)==0 && x==wm.bic,  ri(:,3), 'UniformOutput',0)));
        wm.n_iters=size(ri{wm.ri_num,2},1);
        if isempty(wm.ri_num)==0
            wm.id_num=find(strcmp(rd(:,1), wm.modname).*cell2mat((cellfun(@(x) isempty(x)==0 && size(x,1)==wm.n_iters,  rd(:,3), 'UniformOutput',0))));
            if isempty(wm.id_num)==0
                
                
                % RECORD all details
                nrr(modnum, 1:length(rr(m,:)))=rr(m,:);
                nrd(modnum, :)=rd(wm.id_num,:);
                nri(modnum, :)=ri(wm.ri_num,:);
                modnum=modnum+1; 
                
            else  missingres{length(missingres)+1,1}=wm.modname;
            end
        
        else  missingres{length(missingres)+1,1}=wm.modname;
        end
    else
        missingres{length(missingres)+1,1}=wm.modname;
        
    end
    
    
end


end

