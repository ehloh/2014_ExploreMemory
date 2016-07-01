% Clean up runs (multiple options)
clear all; close all hidden; clc

request.whatsmissing=0;
request.cleanpartialfit= 1;
request.doremove=0;
request.nans_4failedruns=0;
% request.checkparams_inspace=0;


for o=1:1 % General setup 
 % Folders
w=pwd;
if strcmp(w(1), '/')==0;  path(pathdef); addpath('D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling')
else error('Folder!'); end
end

%% % Take a partialfit or a workspace from a crashed fit and clean it up, turn  it into a proper results file
% Removes all fits that (a) outside the specified range, (b) not converged,
% or (c) where r_res, r_iterd and r_iterations don't match up.

if request.cleanpartialfit
    wherefile='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling\1 Inputs\2b Hierar fits\New\v2';
    infilename= 'res_hierarfitmodels_cF (15-Apr-2015)';
    lastokrun=2;
    
    for o=1:1
        ishierar=1;
        
        % ---------------------------------------------------------------------------------------------------
        cd(wherefile);  load([infilename '.mat']); openvar r_res;
        input(['Check r_res:  Last ok run is :' num2str(lastokrun) '  ?   ' ]);
        
        r_res=r_res(1:lastokrun,:);
        for m=1:size(r_res,1); if isempty(r_res{m,2})==1;  r_res{m,3}=nan; end; end % unconverged
        if ishierar ==1;   [ r_res , r_iterations, r_iterd, missingres] = f_hierarmatch(r_res , r_iterations, r_iterd);  openvar missingres, disp('Models that need rerunning are listed in missingres'); end
        %
        for m=1:size(r_res,1)
            ri(m,:)=r_iterations(strcmp(r_iterations(:,1) , r_res{m,1}),:);
            rd(m,:)=r_iterd(strcmp(r_iterd(:,1) , r_res{m,1}),:);
            dm(m,:)= details.models(strcmp(details.models(:, 1), r_res{m,1}), :);
        end
        
        details.models=dm; details.whichmodels=details.models(:,1); details.n_models=length(details.whichmodels);
        r_iterations=ri; if ishierar ==1;   r_iterd=rd; end
        
        for m=1:length(details.whichmodels)  % CHECK
            if  strcmp(details.whichmodels{m,1}, r_res{m,1})+  strcmp(r_res{m,1}, r_iterd{m,1}) +  strcmp(r_res{m,1}, r_iterations{m,1})~=3
                error([r_res{m,1}   '     : Model matchup FAILS']);
            end
        end
        
        % SAVE ---------------------------------------------------------------------------------------------------      
        resfilename=f_newname([infilename ' CleanPartial'], wherefile);
        input(['Save as     ' resfilename  '  ?']);
        if ishierar ==1; save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');   % Save a hierarfit
        else  save(resfilename, 'details', 'r_iterations','r_res', 'errorlog');  % Saved a fixed fit
        end
    end
end

%% Identify models that are still missing
% Checks for what other models need to be run. Assumes list of models here
% is the fulll list.
% Count only models where r_res, r_iterd and r_iterations match up.

if request.whatsmissing
    ishierar=1;
    wherefile='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling\1 Inputs\2b Hierar fits';
    infilename= 'res_hierarfitmodels_cF (13-Apr-2015) v2 bpji7e_L3523';
        
    % ---------------------------------------------------------------------------------------------------
    cd(wherefile);  load([infilename '.mat']); 
    
    for o=1:1 % All models (assume cF)
        
        
        details.ct_modnums=[1 3 4 5 6 11 12 13 17 19 21 23]; % true models for c
        which=1:24;
        
        modfams.b={'b01';'b02_f';'b03_e';'b04_uw';'b05_vw';'b06_ow';'b07_fe';'b08_fuw';'b09_fvw';'b10_fow';'b11_euw';'b12_evw';'b13_eow';'b14_feuw';'b15_fevw';'b16_feow';'b17_yw';'b18_fyw';'b19_eyw';'b20_feyw'; 'b21_kw';'b22_fkw';'b23_ekw';'b24_fekw'};
        modfams.bp={'bp01';'bp02_f';'bp03_e';'bp04_uw';'bp05_vw';'bp06_ow';'bp07_fe';'bp08_fuw';'bp09_fvw';'bp10_fow';'bp11_euw';'bp12_evw';'bp13_eow';'bp14_feuw';'bp15_fevw';'bp16_feow';'bp17_yw';'bp18_fyw';'bp19_eyw';'bp20_feyw'; 'bp21_kw';'bp22_fkw';'bp23_ekw';'bp24_fekw'};
        modfams.bm={'bm01';'bm02_f';'bm03_e';'bm04_uw';'bm05_vw';'bm06_ow';'bm07_fe';'bm08_fuw';'bm09_fvw';'bm10_fow';'bm11_euw';'bm12_evw';'bm13_eow';'bm14_feuw';'bm15_fevw';'bm16_feow';'bm17_yw';'bm18_fyw';'bm19_eyw';'bm20_feyw'; 'bm21_kw';'bm22_fkw';'bm23_ekw';'bm24_fekw'};
        modfams.bpm ={'bpm01';'bpm02_f';'bpm03_e';'bpm04_uw';'bpm05_vw';'bpm06_ow';'bpm07_fe';'bpm08_fuw';'bpm09_fvw';'bpm10_fow';'bpm11_euw';'bpm12_evw';'bpm13_eow';'bpm14_feuw';'bpm15_fevw';'bpm16_feow';'bpm17_yw';'bpm18_fyw';'bpm19_eyw';'bpm20_feyw'; 'bpm21_kw';'bpm22_fkw';'bpm23_ekw';'bpm24_fekw'};
        %
        modfams.bi={'bi01';'bi02_f';'bi03_e';'bi04_uw';'bi05_vw';'bi06_ow';'bi07_fe';'bi08_fuw';'bi09_fvw';'bi10_fow';'bi11_euw';'bi12_evw';'bi13_eow';'bi14_feuw';'bi15_fevw';'bi16_feow';'bi17_yw';'bi18_fyw';'bi19_eyw';'bi20_feyw'; 'bi21_kw';'bi22_fkw';'bi23_ekw';'bi24_fekw'};
        modfams.bpi ={'bpi01';'bpi02_f';'bpi03_e';'bpi04_uw';'bpi05_vw';'bpi06_ow';'bpi07_fe';'bpi08_fuw';'bpi09_fvw';'bpi10_fow';'bpi11_euw';'bpi12_evw';'bpi13_eow';'bpi14_feuw';'bpi15_fevw';'bpi16_feow';'bpi17_yw';'bpi18_fyw';'bpi19_eyw';'bpi20_feyw'; 'bpi21_kw';'bpi22_fkw';'bpi23_ekw';'bpi24_fekw'};
        modfams.bmi ={'bmi01';'bmi02_f';'bmi03_e';'bmi04_uw';'bmi05_vw';'bmi06_ow';'bmi07_fe';'bmi08_fuw';'bmi09_fvw';'bmi10_fow';'bmi11_euw';'bmi12_evw';'bmi13_eow';'bmi14_feuw';'bmi15_fevw';'bmi16_feow';'bmi17_yw';'bmi18_fyw';'bmi19_eyw';'bmi20_feyw'; 'bmi21_kw';'bmi22_fkw';'bmi23_ekw';'bmi24_fekw'};
        modfams.bpmi ={'bpmi01';'bpmi02_f';'bpmi03_e';'bpmi04_uw';'bpmi05_vw';'bpmi06_ow';'bpmi07_fe';'bpmi08_fuw';'bpmi09_fvw';'bpmi10_fow';'bpmi11_euw';'bpmi12_evw';'bpmi13_eow';'bpmi14_feuw';'bpmi15_fevw';'bpmi16_feow';'bpmi17_yw';'bpmi18_fyw';'bpmi19_eyw';'bpmi20_feyw'; 'bpmi21_kw';'bpmi22_fkw';'bpmi23_ekw';'bpmi24_fekw'};
        %
        modfams.bj={'bj01';'bj02_f';'bj03_e';'bj04_uw';'bj05_vw';'bj06_ow';'bj07_fe';'bj08_fuw';'bj09_fvw';'bj10_fow';'bj11_euw';'bj12_evw';'bj13_eow';'bj14_feuw';'bj15_fevw';'bj16_feow';'bj17_yw';'bj18_fyw';'bj19_eyw';'bj20_feyw'; 'bj21_kw';'bj22_fkw';'bj23_ekw';'bj24_fekw'};
        modfams.bpj={'bpj01';'bpj02_f';'bpj03_e';'bpj04_uw';'bpj05_vw';'bpj06_ow';'bpj07_fe';'bpj08_fuw';'bpj09_fvw';'bpj10_fow';'bpj11_euw';'bpj12_evw';'bpj13_eow';'bpj14_feuw';'bpj15_fevw';'bpj16_feow';'bpj17_yw';'bpj18_fyw';'bpj19_eyw';'bpj20_feyw'; 'bpj21_kw';'bpj22_fkw';'bpj23_ekw';'bpj24_fekw'};
        modfams.bjm={'bjm01';'bjm02_f';'bjm03_e';'bjm04_uw';'bjm05_vw';'bjm06_ow';'bjm07_fe';'bjm08_fuw';'bjm09_fvw';'bjm10_fow';'bjm11_euw';'bjm12_evw';'bjm13_eow';'bjm14_feuw';'bjm15_fevw';'bjm16_feow';'bjm17_yw';'bjm18_fyw';'bjm19_eyw';'bjm20_feyw'; 'bjm21_kw';'bjm22_fkw';'bjm23_ekw';'bjm24_fekw'};
        modfams.bpjm={'bpjm01';'bpjm02_f';'bpjm03_e';'bpjm04_uw';'bpjm05_vw';'bpjm06_ow';'bpjm07_fe';'bpjm08_fuw';'bpjm09_fvw';'bpjm10_fow';'bpjm11_euw';'bpjm12_evw';'bpjm13_eow';'bpjm14_feuw';'bpjm15_fevw';'bpjm16_feow';'bpjm17_yw';'bpjm18_fyw';'bpjm19_eyw';'bpjm20_feyw'; 'bpjm21_kw';'bpjm22_fkw';'bpjm23_ekw';'bpjm24_fekw'};
        %
        modfams.bji={'bji01';'bji02_f';'bji03_e';'bji04_uw';'bji05_vw';'bji06_ow';'bji07_fe';'bji08_fuw';'bji09_fvw';'bji10_fow';'bji11_euw';'bji12_evw';'bji13_eow';'bji14_feuw';'bji15_fevw';'bji16_feow';'bji17_yw';'bji18_fyw';'bji19_eyw';'bji20_feyw'; 'bji21_kw';'bji22_fkw';'bji23_ekw';'bji24_fekw'};
        modfams.bpji={'bpji01';'bpji02_f';'bpji03_e';'bpji04_uw';'bpji05_vw';'bpji06_ow';'bpji07_fe';'bpji08_fuw';'bpji09_fvw';'bpji10_fow';'bpji11_euw';'bpji12_evw';'bpji13_eow';'bpji14_feuw';'bpji15_fevw';'bpji16_feow';'bpji17_yw';'bpji18_fyw';'bpji19_eyw';'bpji20_feyw'; 'bpji21_kw';'bpji22_fkw';'bpji23_ekw';'bpji24_fekw'};
        modfams.bjmi={'bjmi01';'bjmi02_f';'bjmi03_e';'bjmi04_uw';'bjmi05_vw';'bjmi06_ow';'bjmi07_fe';'bjmi08_fuw';'bjmi09_fvw';'bjmi10_fow';'bjmi11_euw';'bjmi12_evw';'bjmi13_eow';'bjmi14_feuw';'bjmi15_fevw';'bjmi16_feow';'bjmi17_yw';'bjmi18_fyw';'bjmi19_eyw';'bjmi20_feyw'; 'bjmi21_kw';'bjmi22_fkw';'bjmi23_ekw';'bjmi24_fekw'};
        modfams.bpjmi={'bpjmi01';'bpjmi02_f';'bpjmi03_e';'bpjmi04_uw';'bpjmi05_vw';'bpjmi06_ow';'bpjmi07_fe';'bpjmi08_fuw';'bpjmi09_fvw';'bpjmi10_fow';'bpjmi11_euw';'bpjmi12_evw';'bpjmi13_eow';'bpjmi14_feuw';'bpjmi15_fevw';'bpjmi16_feow';'bpjmi17_yw';'bpjmi18_fyw';'bpjmi19_eyw';'bpjmi20_feyw'; 'bpjmi21_kw';'bpjmi22_fkw';'bpjmi23_ekw';'bpjmi24_fekw'};
        %
        allmods=[
            modfams.b(which)
            modfams.bp(which)
            modfams.bm(which)
            modfams.bpm(which)
            modfams.bi(which)   %
            modfams.bpi(which)
            modfams.bmi(which)
            modfams.bpmi(which)
            modfams.bj(which)   %
            modfams.bpj(which)
            modfams.bjm(which)
            modfams.bpjm(which)
            modfams.bji(which)  %
            modfams.bpji(which)
            modfams.bjmi(which)
            modfams.bpjmi(which)
            ];
    end
    
    % Check for models that are missing modelspaces
    if ishierar ==1;   [ r_res , r_iterations, r_iterd, missingres] = f_hierarmatch(r_res , r_iterations, r_iterd);  openvar missingres, disp('Models that need rerunning are listed in missingres'); end
    
    %
    r_modthere=zeros(length(allmods),1);
    for m=1:length(allmods)
        if sum(strcmp(r_res(:,1), allmods{m}))==1 && isempty(r_res(strcmp(r_res(:,1), allmods{m}),2))==0
            r_modthere(m)=1;
        else r_modthere(m)=0;
        end
    end
    if sum(r_modthere)==length(allmods); disp('All mods there! :)'); 
    else missingmods=allmods(r_modthere==0); openvar missingmods, disp('Missing mods:'), disp(missingmods)
    end
    
    
    
end


%% Remove models with certain parameters from the model space
% OR, remove specific models
% ALSO, remove duplicates (leave removepars & removemods empty)


if request.doremove
    removepars={'i'};
    removemods={};
    
    wherefile='D:\Dropbox\SANDISK\8 Explore Mem\3b Modelling\1 Inputs\2b Hierar fits\New\v1';
    filename='res_hierarfitmodels_cF (15-Apr-2015) Diego1';
    ishierar=1;
    for o=1:1
        % ------------------------------------------------------
        
        cd(wherefile);    load([filename '.mat']);
        
        for p=1:length(removepars)
            
            for m=1:size(r_res,1)  % res
                if isempty(strfind(r_res{m,1}, removepars{p}))==0
                    r_res{m,1}='zzzzzzzzzzzz';
                end
            end
            for m=1:size(r_iterations,1) % r_iterations
                if isempty(strfind(r_iterations{m,1}, removepars{p}))==0
                    r_iterations{m,1}='zzzzzzzzzzzz';
                end
            end
            if ishierar ==1;
                for m=1:size(r_iterd,1) % r_iterd
                    if isempty(strfind(r_iterd{m,1}, removepars{p}))==0
                        r_iterd{m,1}='zzzzzzzzzzzz';
                    end
                end
            end
            
            % DETAILS #############
            for m=1:size(details.models,1)   % details.models
                if isempty(strfind(details.models{m,1}, removepars{p}))==0
                    details.models{m,1}='zzzzzzzzzzzz';
                end
            end
            
            
        end
        for o=1:1  % Remove duplicates
            newfits={}; k=1;
            for m=1:size(r_res ,1)  % res
                if sum(strcmp(r_res(:,1), r_res{m,1}))>1 & strcmp( r_res{m,1},  'zzzzzzzzzzzz')==0
                    wr.fits=r_res(find(strcmp(r_res(:,1), r_res{m,1})),:);
                    %             wr.fits=
                    r_res(find(strcmp(r_res(:,1), r_res{m,1})),1) = repmat(    {'zzzzzzzzzzzz'}, sum(strcmp(r_res(:,1), r_res{m,1})), 1);
                    
                    %
                    wr.Ls=cell2mat(wr.fits(:, 3));
                    if sum(wr.Ls==min(wr.Ls))~=1; error('MULTIPLE matches. Just delete one manually'); end
                    wr.bestfit=wr.fits ( find(wr.Ls==min(wr.Ls)  ), :);
                    newfits(k,  1:  size(wr.fits,2)) =  wr.bestfit;   k=k+1;
                end
            end
            r_res=[r_res; newfits];
            
            
            
            for m=1:size(details.models,1)  % res
                if sum(strcmp(details.models(:,1), details.models{m,1}))>1 & strcmp( details.models{m,1},  'zzzzzzzzzzzz')==0
                    details.models{m,1}='zzzzzzzzzzzz';
                end
                
            end
            
        end
        
        
        % TO DO: Mark models to remove!!
        
        % Get rid of
        r_res=sortrows(r_res,1); r_res(strcmp(r_res(:,1), 'zzzzzzzzzzzz'), :)=[];
        r_iterations=sortrows(r_iterations,1);  r_iterations(strcmp(r_iterations(:,1), 'zzzzzzzzzzzz'), :)=[];
        if ishierar ==1;   r_iterd= sortrows(r_iterd,1);  r_iterd(strcmp(r_iterd(:,1), 'zzzzzzzzzzzz'), :)=[]; end
        details.models=sortrows(details.models,1); details.models(strcmp(details.models(:,1), 'zzzzzzzzzzzz'), :)=[];
        details.whichmodels=sortrows(details.models(:,1),1);    details.n_models=length(details.whichmodels);
        %     r_res=sortrows(r_res,3);
        
        % SAVE ---------------------------------------------------------------------------------------------------
        if ishierar ==1  ; resfilename=['FIXED res_hierarfitmodels_' details.task ' (' date ')']; % Save a hierarfit
        else  resfilename=['FIXED res_fitmodels_' details.task ' (' date ')'];    % Saved a fixed fit
        end
        resfilename=f_newname(resfilename, wherefile);
        input(['Save as     ' resfilename  '  ?']);
        if ishierar ==1; save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');   % Save a hierarfit
        else  save(resfilename, 'details', 'r_iterations','r_res', 'rc','errorlog');  % Saved a fixed fit
        end
    end
end

%% Fill in NANs for failed runs, identify failed runs (sent to back), sort by best run

if request.nans_4failedruns
    wherefile='/Users/EleanorL/Dropbox/sandisk/8 explore mem/3b modelling/1 Inputs/2b Hierar fits/Hierarchical exhaustive both';
    filename='res_hierarfitmodels_cF (13-Apr-2015)1';
    ishierar =1;
    
    for o=1:1
        %     --------------------------------------------------------------------------------------------------------
        fillwhat=nan;
        %
        cd(wherefile), load([filename '.mat']);
        for m=1:size(r_res,1) % Mark bad mods
            if isempty(r_res{m,3}) || r_res{m,3}==999 || r_res{m,3}<0
                r_res{m,3}= fillwhat;
            end
        end
        
        % Sort
        r_res=sortrows(r_res,3);
        badmods=r_res(find(isnan(cell2mat(r_res(:, 3)))),1);
        
        
        % SAVE ---------------------------------------------------------------------------------------------------
        if ishierar ==1  ; resfilename=['FLAGNonConverg res_hierarfitmodels_' details.task ' (' date ')']; % Save a hierarfit
        else  resfilename=['FLAGNonConverg res_fitmodels_' details.task ' (' date ')'];    % Saved a fixed fit
        end
        resfilename=f_newname(resfilename, wherefile);
        input(['Save as     ' resfilename  '  ?']);
        if ishierar ==1; save(resfilename, 'details', 'r_iterations','r_iterd', 'r_res', 'rc','errorlog');   % Save a hierarfit
        else  save(resfilename, 'details', 'r_iterations','r_res', 'rc','errorlog');  % Saved a fixed fit
        end
        
        
        % Display + next steps
        disp('BAD models that did not converge:'), disp(badmods), openvar badmods
        %     edit a3_hierarfit_models.m  % PUT THEM HERE AND RUN THEM!
    end
end

