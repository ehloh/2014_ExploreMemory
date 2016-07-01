% Check for matches between rr_res, rr_iterd, rr_iterations
% MAke sure that that fit that is stored in r_res corresponds to r_iterd
% and r_iterations. Is important for resamping, etc/ 
clear all; close all hidden; clc


wheredata='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling\1 Inputs\2b Hierar fits';
filename='res_hierarfitmodels_cF (13-Apr-2015) v1 bpji16_L3090';
% filename='res_hierarfitmodels_cF (13-Apr-2015) v2 bpji7e_L3523';
% filename='res_hierarfitmodels_cF (13-Apr-2015) both bpji14u_L6618';



%% Remove non-converged and those that are in need of patching 

cd(wheredata), load([filename '.mat']);

d_noconverge={}; d_needpatch={};  % Name, rr_res details, rr_iterations details, rr_iterd, MESSAGE
rr_res=sortrows(r_res,1); rr_iterd=sortrows(r_iterd,1); rr_iterations=sortrows(r_iterations,1); details.models=sortrows(details.models,1);
r_res={}; r_iterd={}; r_iterations={}; mnum=1;  dm={}; 
for m=1:size(rr_res,1)  % m=res, mi=rr_iterations, md=rr_iterd
    if isempty(rr_res{m,2})==0
        wm.iterations_num=find(cell2mat(cellfun(@(x) isempty(x)==0 && x==rr_res{m,3},rr_iterations(:,3),'UniformOutput',0)));
        if length(wm.iterations_num)==1
            wm.n_iter=size(rr_iterations{wm.iterations_num,2},1);            
            wm.iterd_num=find(strcmp(rr_iterd(:,1),rr_res{m,1}).*cellfun(@(x)size(x,1), rr_iterd(:,3))==wm.n_iter);
            
            
            % RECORD ALL cORRECT DETAILS
            r_res(mnum,:)=rr_res(m,:);           
            r_iterd(mnum,:)=rr_iterd(wm.iterd_num,:);
            r_iterations(mnum,:)=rr_iterations(wm.iterations_num,:);
            dm(mnum,:)=details.models(strcmp( details.models(:,1), rr_res{m,1}),:);
            
            
            mnum=mnum+1; 
        else
            d_needpatch{size(d_needpatch,1)+1,1}=rr_res{m,1};
            d_needpatch{end,2}=rr_res(m,:);
        end
    else
        d_noconverge{length(d_noconverge)+1,1}=rr_res{m,1};
    end
end

% Check
for m=1:size(r_res,1)
    if r_res{m,3} == r_iterations{m,3} && size(r_iterations{m,2},1)==size(r_iterd{m,3},1)
        disp([r_res{m,1} '  ok'])
    else disp([r_res{m,1} '  ARGH'])
    end
end


details.models=dm;
details.whichmodels=details.models(:,1);
details.n_models=length(details.whichmodels);


openvar d_noconverge, openvar  d_needpatch

%
input('SAVE?');

resfilename=[filename ' PATCHok'];
resfilename=f_newname(resfilename, pwd);
save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog'); 

