% Setup data for modelling
clear all; close all hidden; clc

% Do what?
details.v1_1_ok={'p01_AF' ; 'p01_AS' ; 'p02_TY' ; 'p03_JV' ; 'p06_MK' ; 'p08_AK' ; 'p10_SG' ; 'p11_OG' ; 'p30_AR' ; 'p31_RD' ; 'p33_AA' ; 'p34_BE'};
details.v1_2_ok={'p14_CP';'p15_SZ';'p23_SG';'p24_VP';'p26_PB';'p27_PV';'p29_MY';'p35_EO';'p37_ND';'p38_SL';'p39_JH';'p40_CC'; 'p45_SO';};
details.taskversion='Both';
details.specificsubjects=[details.v1_1_ok; details.v1_2_ok];

for o1=1:1 % General setup
%     where.root='D:\Dropbox\SANDISK\';
    where.root='/Users/EleanorL/Dropbox/SANDISK/';
    where.data=[where.root '8 Explore Mem' filesep '2 All data' filesep 'Both'];
    where.analysisscript=[where.root filesep '8 Explore Mem' filesep '3a Analysis scripts'];
    path(pathdef); addpath([where.root '8 Explore Mem']); addpath(where.analysisscript);
    
    % Fetch subjects
    [w.s w.s1 details.koshertable]=xlsread([where.analysisscript filesep 'i_subjectsok.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [details.subjects details.n_subjs details.datalog] = f_selectsubjects(details.koshertable, details.specificsubjects, details.koshertable, 'Include');

    disp(' ##################################################')
    disp(['No. of subjects:  ' num2str(details.n_subjs)])
    if isempty(details.specificsubjects)==0
        disp(details.subjects)
    end
    input('Hit enter to continue     '); disp(' ')
    disp(' ##################################################')
end

%% Load subject data + stick into column format for modelling 
%   Variable 'subjdata': Col 1=Subject, Col 2=Conflict task

for o1=1:1 % Columns in modelling data
    
    % Columns for modelling data
    col.Choice=3;
    col.EnvThreat=6;
    col.NTokens=2;
    col.pLoss=1;
    col.Entropy=4;
    col.VExplore=5;
    col.EV=10;
    col.OutcomeMagnitude=11;
    col.OutcomeMean=12;
    col.OutcomeVariance=13;
    col.Task=9;
    col.Trialnum=8;
    col.EntropyNTok=7;
%     col.EntropyEV=14;
    col.StanDev=16;
%     col.vMeanVar=17;
    
    % 
    details.col=col;
end


% Load data
subjdata=cell(details.n_subjs, 4);
for s=1:details.n_subjs
    
    ws=load([where.data filesep details.subjects{s} filesep details.subjects{s} '_file_3ExploreEnc.mat']);
    if s==1; col.f=ws.cF.col; end
    
    ws.ts=ws.cF.data(ws.cF.data(:,col.f.Trialvalid)==1, :);
    ws.d=nan*zeros(size(ws.ts,1), 7);
    ws.d(:,col.NTokens)= ws.ts(:, col.f.NTokenPairs) *2;
    ws.d(:,col.EnvThreat)=ws.ts(:, col.f.EnvThreat)/4;
    ws.d(:,col.Choice)=ws.ts(:, col.f.Resp1);
    ws.d(:,col.Trialnum)=  1:size(ws.ts,1);
    ws.d(:,col.Task)=1;
    ws.d(:,col.OutcomeMagnitude)=ws.ts(:,col.f.OutcomeMag);
    
    % Write parameters
    ws.cf=ws.d(ws.d(:,col.Task)==1, :);
    [ ws.cf] = fpar_conflict( ws.cf, col);
    
    for e=1:4
        for n=1:4
            wc=ws.cf(ws.cf(:, col.EnvThreat)*4==e & ws.cf(:, col.NTokens)/2==n, :);
            for c=1:3
                mplot{s,c}(5-e,n)=mean(wc(:, col.Choice)==c);
                
            end
        end
    end
    
    
    subjdata{s, 1}=details.subjects{s};
    subjdata{s, 2}=ws.cf;
    subjdata{s, 3}=[];
    subjdata{s, 4}=ws.cf;
    
    ws=[];
end


% Check plots?
% s=3;
% c=1; subplot(1,3,c); imagesc(mplot{s,c}), colorbar, axis square, caxis([0 1])
% c=2; subplot(1,3,c); imagesc(mplot{s,c}), colorbar, axis square, caxis([0 1])
% c=3; subplot(1,3,c); imagesc(mplot{s,c}), colorbar, axis square, caxis([0 1])
% 
% 
% error


%%

input('Save?');
save([pwd  filesep '1 Inputs' filesep 'All data (' date ').mat'], 'subjdata', 'details')



