% Establish subject-specific parameters to be run in other scripts
clear all; close all; clc

testing=0;
for o1=1:1 % Documentation
    % This script sets up columns 1-6; executed in other scripts during tasks themselves
    %         Col 1:       TrialType
    %         Col 2:       EnvThreat (1=Low, 2=High)
    %         Col 3:       NToken pairs (1 or 3)
    %         Col 4:       Bomb present?
    %         Col 5:       Activated bomb present?
    %         Col 6:       Trial no.
    %         Col 7:       Monitoring trial? (Learning phase only)
    %         -
    %         Col 8:       [1st response]  Accept, Reject or Explore (1=Accept, 2=Reject, 3=Explore)
    %         Col 9:       [1st response]  RT
    %         Col 10:      [2nd response]  Accept or Reject (1=Accept, 2=Reject, 3=Explore)
    %         Col 11:      [2nd response]  RT
    %         Col 12:      [Exploration] Did exploration reveal a bomb?  (1=Yes, 0=No, 999=Not explored)
    %         Col 13:      Valid trial?
    %         Col 14:      Outcome magnitude (NTokens)
    %         Col 15:     Item stim no.
    %         Col 16:     Item type (1=Manmade, 2=Natural)
    
    col.Trialtype=1;
    col.EnvThreat=2;
    col.NTokenPairs=3;
    col.BombPresent=4;
    col.BombActivated=5;
    col.Trialnum=6;
    col.MonitorTask=7;
    col.Resp1=8;
    col.RT1=9;
    col.Resp2=10;
    col.RT2=11;
    col.BombExplored=12;
    col.Trialvalid=13;
    col.OutcomeMag=14;
    col.ItemStim=15;
    col.ItemType=16;
    col.ItemPos=17;
end
for o1=1:1 % Execution
    
    if testing==1
        p.subject=input('Subject ID:   ', 's');
    else
        p.subject='t1';
    end
end
for o1=1:1 % General settings
    
    % Design (Encoding phase)
    p.design.nreps_percell=16;
    p.design.EnvThreat_Levels=4;
    p.design.NTokenPairs_Levels=4;
    p.design.NSpaces_TokenPairs=4;
    p.nTrials=p.design.nreps_percell*p.design.EnvThreat_Levels*p.design.NTokenPairs_Levels;
    p.design.ExplorationCost=2; % In NTokens, not pairs
    p.design.FixedLoss=-8;
    %
    p.nItemsEncoded=p.nTrials;
    p.nItemsFoil=ceil(p.nItemsEncoded/2);
    p.nItemsAll=p.nItemsEncoded+p.nItemsFoil;
    
    % Env-Learning specs
    p.learn.nreps_percell=10;
    p.learn.EnvThreat_Levels=p.design.EnvThreat_Levels;
    p.learn.NTokenPairs_Levels=p.design.NTokenPairs_Levels;
    p.learn.NSpaces_TokenPairs=p.design.NSpaces_TokenPairs;
    p.learn.MonitorTaskPercent=0.1;
    
    
    % Display
    p.disp.tokenY_row1=50; % Display
    p.disp.tokenY_row2=-50;
    p.disp.tokenpos=[-150; -50; 50; 150];
    p.disp.tokensize=80;
    p.disp.tokengap=60;
    p.disp.settokencolor=[0.8 0.8 0.8];    
    p.disp.set_contingencycolours={[0.1 0.3 0.7] 'Blue'; [0.1 0.4 0.4] 'Green'; [0.7 0.4 0] 'Orange';[0.4 0 0.6] 'Purple';[0.7 0.7 0.0] 'Yellow';[0.8 0.5 0.5] 'Pink'};
    p.disp.color_bomb=[0.7 0 0];
    p.disp.color_nobomb=[0 0.7 0];
    p.disp.color_nowinnolose=[1 1 0]; 
    %
    p.disp.ItemSize=250;
    p.disp.ItemPos1_x=-180;
    p.disp.ItemPos1_y=160;
    p.disp.ItemPos2_x=180;
    p.disp.ItemPos2_y=160;
    p.disp.ItemPos3_x=-180;
    p.disp.ItemPos3_y=-80;
    p.disp.ItemPos4_x=180;
    p.disp.ItemPos4_y=-80;
        
    % Timings
    p.disp.time_firstofferonly=2000;
    p.disp.time_firstofferitem=2000;
    p.disp.time_secondoffer=p.disp.time_firstofferonly;
    p.disp.time_outcome=1000;
    
    % Keypresses
    p.keys.left1=26; % z
    p.keys.left2=24; % x
    p.keys.left3=3; % c
    p.keys.right1=97; 
    p.keys.right2=100;
    p.keys.right3=98;
    
    % Misc
    rand('state',sum(100*clock));
end

% Generate Subject-specific parameters
p.disp.set_contingencycolours(:,3)=num2cell(rand(size(p.disp.set_contingencycolours,1),1));
p.disp.contingencycolours=sortrows(p.disp.set_contingencycolours,3);

% Env-Learning phase parameters
[ learnpar learnnorm p.learn] = f_generate_taskstruc(p.learn, col);
learnpar(:,col.MonitorTask)=zeros(size(learnpar,1),1);
learnpar(:,col.Trialnum)=rand(size(learnpar,1),1); learnpar=sortrows(learnpar, col.Trialnum);
learnpar(1:floor(p.learn.MonitorTaskPercent*size(learnpar,1)),col.MonitorTask)=1;
learnpar(:,col.Resp1)=1;learnpar(:,col.Resp2)=nan;
for e=1:p.design.EnvThreat_Levels
    for n=1:p.design.NTokenPairs_Levels
        learnpar(learnpar(:,col.EnvThreat)==e & learnpar(:,col.NTokenPairs)==n, col.Trialtype)=(e-1)*p.design.NTokenPairs_Levels+n;
    end
end

% Encoding phase parameters
[ trialpar norm p.design] = f_generate_taskstruc(p.design, col);
trialpar(trialpar(:,col.EnvThreat)==1 & trialpar(:,col.NTokenPairs)==p.design.NTokenPairs_Npairs(1), col.Trialtype)=1;
trialpar(trialpar(:,col.EnvThreat)==1 & trialpar(:,col.NTokenPairs)==p.design.NTokenPairs_Npairs(2), col.Trialtype)=2;
trialpar(trialpar(:,col.EnvThreat)==2 & trialpar(:,col.NTokenPairs)==p.design.NTokenPairs_Npairs(1), col.Trialtype)=3;
trialpar(trialpar(:,col.EnvThreat)==2 & trialpar(:,col.NTokenPairs)==p.design.NTokenPairs_Npairs(2), col.Trialtype)=4;
for e=1:p.design.EnvThreat_Levels
    for n=1:p.design.NTokenPairs_Levels
        trialpar(trialpar(:,col.EnvThreat)==e & trialpar(:,col.NTokenPairs)==n, col.Trialtype)=(e-1)*p.design.NTokenPairs_Levels+n;
    end
end

% Add item stimuli
load(['Stimuli' filesep 'stimlist.mat']);
itemlist(:,4)=num2cell(rand(size(itemlist,1),1)); % Randomize assignment to encoded or foil (control for item type)
itemlist=sortrows(itemlist,4);
wi.t1=itemlist(cell2mat(itemlist(:,2))==1,:);
wi.t2=itemlist(cell2mat(itemlist(:,2))==2,:);
wi.t1=wi.t1(1:p.nItemsAll/2,:);
wi.t2=wi.t2(1:p.nItemsAll/2,:);
wi.enc=[cell2mat([wi.t1(1:p.nItemsEncoded/2,[3 2]); wi.t2(1:p.nItemsEncoded/2,[3 2])]) rand(p.nItemsEncoded,1)];
wi.enc=sortrows(wi.enc,3);
p.ItemFoils=[wi.t1(p.nItemsEncoded/2+1:end,2:3); wi.t2(p.nItemsEncoded/2+1:end,2:3)];
p.ItemEnc=wi.enc(:,[2 1]);
trialpar(:, [col.ItemStim col.ItemType])=wi.enc(:, 1:2);
trialpar(:,col.ItemPos)=randi(4,p.nTrials,1);

% Randomization & setup into tasks
learnpar(:,col.Trialnum)=rand(size(learnpar,1),1);  % Randomize LearnEnv
learnpar=sortrows(learnpar, col.Trialnum); learnpar(:,col.Trialnum)=1:size(learnpar,1);
data_learnenv=nan*zeros(size(learnpar,1),14);
data_learnenv(:,1:size(learnpar,2))=learnpar;
trialpar(:,col.Trialnum)=rand(p.nTrials,1); % Randomize conflict
trialpar=sortrows(trialpar, col.Trialnum); trialpar(:,col.Trialnum)=1:size(trialpar,1);
data_cF=nan*zeros(size(trialpar,1),14);
data_cF(:,1:size(trialpar,2))=trialpar;
% trialpar(:,col.Trialnum)=rand(p.nTrials,1); % Randomize control
% trialpar=sortrows(trialpar, col.Trialnum); trialpar(:,col.Trialnum)=1:size(trialpar,1);
% data_ct=nan*zeros(size(trialpar,1),14);
% data_ct(:,1:size(trialpar,2))=trialpar;

% Save 
Learn.par=p;
Learn.col=col;
Learn.learnpar=data_learnenv;
cF.par=p;
cF.col=col;
cF.trialpar=data_cF;
% ct.par=p;
% ct.col=col;
% ct.trialpar=data_ct;
save(['Data' filesep p.subject '_file_1taskpar.mat'], 'cF', 'Learn');
try % Transfer file to DeletedDaily 
    movetowhere='\\Asia\DeletedDaily\EL ExploreMem Data';
    if isdir(movetowhere)==0; mkdir(movetowhere);end
    copyfile(['Data' filesep p.subject '_file_1taskpar.mat'],  [movetowhere filesep p.subject '_file_1taskpar.mat']);
    w.warning='Parameters created successfully';
catch
    w.warning='WARNING: Data not transfered to Deleted Daily. Transfer files manually';
end
disp('Parameters created successfully'); disp(' ')

% Display to experimenter
disp('----------------------------------------------')
disp('NOTICE: Subject turn away!')
input('Experimenter hit enter to proceed   '); disp(' ')
disp('Colours in ascending order of threat:'); disp(' ')
disp(p.disp.contingencycolours(:,2));
disp('----------------------------------------------')


%% Generate plots


% To plot the design!
plotdesign=1;


if plotdesign
    disp('Go to specified path, find fpar_conflict, and turn on the plots');
    addpath('D:\Dropbox\SANDISK\8 Explore Mem\3 Analysis scripts');
    col.NTokens=col.NTokenPairs;
    col.pLoss=18;
    col.Entropy=19;
    col.EV=20;
    col.EntropyNTok=21;
    
    [ dplot] = fpar_conflict(trialpar,  col);
end