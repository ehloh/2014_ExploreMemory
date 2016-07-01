% Collate misc data, save to subjects' Encoding & Memory files (for printout to data tables)
clear all; close all hidden; clc

where.where='F:\8 Explore Mem';
% where.where='/Volumes/SANDISK/8 Explore Mem';

%
log.TaskVersion='v1_1';


for o1=1:1 % General setup  
    
    where.alldata=[where.where filesep '2 All Data'];
    where.execscript=([where.where filesep '1 Explore Mem Task']); addpath(where.execscript)
    addpath([where.where filesep '3 Analysis scripts'])
    
    % Fetch subjects
    log.specificsubjects=[];
    [w.s w.s1 log.koshertable]=xlsread([where.where filesep '3 Analysis scripts' filesep 'i_subjectsok.xlsx']); % If script cannot find specific subjects here, these subjects need to be excluded on excel sheet
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.koshertable, log.specificsubjects, log.koshertable, 'Include');
       
end

%% Load all subject data (and append details)

disp('Loading subject data ----------------------------------')
subjdata=cell(log.n_subjs,3); % Col 1= Name, Col 2= Trialstats, Col 3=Choice grids, Row n+1 = Group
for s=1:log.n_subjs
    disp(['Subject ' num2str(s) ' (' log.subjects{s} ') -- '])
    ws.e=load([where.alldata filesep log.TaskVersion filesep log.subjects{s} filesep log.subjects{s} '_file_3ExploreEnc.mat']);
    ws.m=load([where.alldata filesep log.TaskVersion filesep log.subjects{s} filesep log.subjects{s} '_file_4Memtest.mat']);
    % ------------------------------------------------------------------------------------------------------------------------

    % (1) TPQ
    
    % when abigail gives me the data, load it, append to 
    
    
    
    
    
    
    % ------------------------------------------------------------------------------------------------------------------------
    % Copy to memstats, save both
%     ws.m.mem.behstats=
    
    
    %
    ws=[];
end

