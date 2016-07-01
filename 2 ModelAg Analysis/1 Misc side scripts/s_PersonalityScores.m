% Read personality scores from excel
clear all; clc


where.data='/Users/EleanorL/Desktop/Participant personality tests';


log.subjects={'p02_TY';'p03_JV'}; % Subject labels (including excel names) must be PRECISE!
log.n_subjs=length(log.subjects);


%% TPQ

% Abigail: Excels need to be saved in xlsx format, then run your eye down
% and make sure that responses are in capital Ts and Fs. If some subjects
% use other formats, may need to adapt script - or do all subjects who use
% one format first, then do other subjects after adapting the script (can specify which subjects
% specifically above). Or, re-enter subjects responses in the necessary
% format (but be careful!!!). 
% Excel file names also need to be specified very precisely - eg P01 rather
% than p01 will cause an error! Can also delete all the subject names from
% the excel filenames, and change the code below. 

% Put TPQ into a format that the TPQ script can then read
for s=1:log.n_subjs
    [ws.n ws.t ws.raw]= xlsread([where.data filesep log.subjects{s} filesep 'TPQ_Questionnaire .xlsx']);
    ws.resp=ws.raw(:,4); % raw responses from excel
    ws.answers=zeros(100,4);
    
    q=1; % Assumes questions run sequentially and no other responses are made
    for i=1:length(ws.resp)
        if iscellstr(ws.resp(i))==1 && strcmp(char(ws.resp(i)), 'T')==1 % || iscellstr(ws.resp(i))==1 && strcmp(char(ws.resp(i)), 'TRUE')==1
            ws.answers(q,1)=q;
            ws.answers(q,4)=1;
            q=q+1;
        elseif iscellstr(ws.resp(i))==1 && strcmp(char(ws.resp(i)), 'F')==1
            ws.answers(q,1)=q;
            ws.answers(q,4)=2;
            q=q+1;
        end
    end
    
    % Checks
    if q>101;
        error(['Check Excel sheet for subjects ' log.subjects{s} '. Too many responses - only one response, in correct format, per question!')
    elseif q<101
        error(['Check Excel sheet for subjects ' log.subjects{s} '. Missing responses - must have at least one response, in correct format, per question!')
    end
    
    % Save in format
    TPQdata=[[0 0 0 0]; ws.answers];
    save([where.data filesep log.subjects{s} filesep  log.subjects{s} '_TPQdata.mat'], 'TPQdata');
    ws=[]; TPQdata=[];
end

%% ImpSS_short

for s=1:log.n_subjs
    [ws.n ws.t ws.raw]= xlsread([where.data filesep log.subjects{s} filesep 'ImpSS_short.xls'], 'Scoring');
    ws.scores=ws.raw([1 2 3 4 8 9 10 11 12 13], :);

    % Append to encoding stage variable
    ws.e=load([where.data filesep log.subjects{s} filesep log.subjects{s} '_file_3ExploreEnc.mat']);
    
    for i=1:size(ws.scores,1)
        eval(['ws.cF.behstats.'  ws.scores{i,1} '=' ws.scores{i,2} ';'])
    end
        
    % Save!!
    cF=[]; cF=ws.cF;
    
    % save
    
    ws=[];
end
