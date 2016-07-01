% Memory test (trial parameters set up 1st script)
%   To alter trial parameters, do so in 1st script
clear all; close all hidden; clc

testing=0;

for o1=1:1 % Documentation (memtest data)
    
    col.EnvThreat=1;
    col.NTokenPairs=2;
    col.BombActivated=3;
    col.Choice=4;
    col.BombExplored=5;
    col.OutcomeMag=6;
    col.Enc_Trialnum=7;
    col.Enc_TypeOK=8;
    %
    col.Item_OldNew=9;           % 1=Old, 2=New
    col.ItemStim=10;
    col.ItemType=11;        % 1=Man-made, 2=Natural
    col.ItemPos=12;
    %
    col.Mem_Trialnum=13;
    col.OldNew=14;          % 1= Old, 2= New
    col.RemKnow=15;         % 1= Remember, 2= Know
    col.SureGuess=16;       % 1= Sure, 2= Guess
    col.Pos=17;
    col.Roc=18;                 % 1=SureNew, 2=GuessNew, 3=GuessKnow, 4=SureKnow, 5=GuessRem, 6=SureRem
    col.CorrectRecognition=19;
    col.CorrectPosition=20;
end
for o1=1:1 % Execution
    
    if testing==1
        p.subject=input('Subject ID:   ', 's');
        w.fullscreen=1;
    else
        p.subject='t1';
        w.fullscreen=0;
    end
end
for o1=1:1 % General settings + load subject's trial parameters
    pp=load(['Data' filesep p.subject '_file_3ExploreEnc.mat']);
    p.disp.ItemSize=pp.cF.par.disp.ItemSize;    
    p.keys.left=97;
    p.keys.right=98;
    p.keys.k1=91;
    p.keys.k2=97;
    p.keys.k3=100;
    p.keys.k4=98;
    
    % Assemble Trial stats
    encitems=nan*zeros(size(pp.cF.rdata,1), 17);
    encitems(:,col.Item_OldNew)=1;
    encitems(:,col.EnvThreat)=pp.cF.rdata(:,pp.cF.col.EnvThreat);
    encitems(:,col.NTokenPairs)=pp.cF.rdata(:,pp.cF.col.NTokenPairs);
    encitems(:,col.BombActivated)=pp.cF.rdata(:,pp.cF.col.BombActivated);
    encitems(:,col.Choice)=pp.cF.rdata(:,pp.cF.col.Resp1);
    encitems(:,col.BombExplored)=pp.cF.rdata(:,pp.cF.col.BombExplored);
    encitems(:,col.OutcomeMag)=pp.cF.rdata(:,pp.cF.col.OutcomeMag);
    encitems(:,col.Enc_Trialnum)=pp.cF.rdata(:,pp.cF.col.Trialnum);
    encitems(:,col.Enc_TypeOK)=pp.cF.rdata(:,pp.cF.col.Trialvalid);
    encitems(:,col.ItemStim)=pp.cF.rdata(:,pp.cF.col.ItemStim);
    encitems(:,col.ItemType)=pp.cF.rdata(:,pp.cF.col.ItemType);
    encitems(:,col.ItemPos)=pp.cF.rdata(:,pp.cF.col.ItemPos);
    %
    foils=nan*zeros(size(pp.cF.par.ItemFoils,1), 17); % Assemble foils
    foils(:,col.Item_OldNew)=2;
    foils(:,col.ItemStim)=cell2mat(pp.cF.par.ItemFoils(:,2));
    data=vertcat(encitems,foils);
    p.nTrials=size(data,1);
    data(:,col.Mem_Trialnum)=rand(p.nTrials,1);
    data=sortrows(data, col.Mem_Trialnum); data(:,col.Mem_Trialnum)=1:p.nTrials;
    
    % Items
    load(['Stimuli' filesep 'stimlist.mat']); p.itemlist=itemlist;
    
    % Misc
    w.line=50;
    w.back=[0.8 0.8 0.8];
    w.res=2;
    rand('state',sum(100*clock));
    
end

% COGENT
config_display(w.fullscreen,w.res, [0 0 0], [1 1 1], 'Helvetica', 40, 9,0, 0); % marker 5g start
config_keyboard;
start_cogent % if using testing laptops, w.res in config dis must be 3! also, w.lineoptions =-330;
cgtext('The computer will now give you', 0, w.line*4)
cgtext('instructions for this task. ', 0, w.line*3)
cgtext(' If the trials start without you getting', 0, w.line*1)
cgtext('instructions, please tell the experimenter', 0, w.line*0)
cgtext(' immediately', 0, w.line*-1)
cgtext('Press any key to continue with the instructions', 0, w.line*-3)
t.startcogent=cgflip(0.8,0.8,0.8);
p.log.start=clock;
waitkeydown(inf);

%% PREP

for o1=1:1 % Set up stock displays 
    
    % Sprite allocations
    sprite.Resp_OldNew=5;
    sprite.Resp_RemKnow=6;
    sprite.Resp_SureGuess=7;
    sprite.Resp_ItemPos=8;
    sprite.ItemPos=9;
    cgpencol([0 0 0]);  cgfont('Helvetica',30)
    
    % Response options
    cgmakesprite(sprite.Resp_OldNew,700,100, [1 0 0]) % OLD/NEW
    cgtrncol(sprite.Resp_OldNew,'r')
    cgsetsprite(sprite.Resp_OldNew)
    cgrect(-200,0,170, 60, [0.3 0.3 0.3])
    cgrect(200,0,170, 60, [0.3 0.3 0.3])
    cgtext('OLD',-200,0)
    cgtext('NEW',200,0)
    cgsetsprite(0)
   
    cgmakesprite(sprite.Resp_RemKnow,700,100, [1 0 0 ]) % REM/KNOW
    cgtrncol(sprite.Resp_RemKnow,'r')
    cgsetsprite(sprite.Resp_RemKnow)
    cgrect(-200,0,170, 60, [0.3 0.3 0.3])
    cgrect(200,0,170, 60, [0.3 0.3 0.3])
    cgtext('REMEMBER',-200,0)
    cgtext('KNOW',200,0)
    cgsetsprite(0)
    cgmakesprite(sprite.Resp_SureGuess,700,100, [1 0 0]) % SURE/GUESS
    cgtrncol(sprite.Resp_SureGuess,'r')
    cgsetsprite(sprite.Resp_SureGuess)
    cgrect(-200,0,170, 60, [0.3 0.3 0.3])
    cgrect(200,0,170, 60, [0.3 0.3 0.3])
    cgtext('SURE',-200,0)
    cgtext('GUESS',200,0)
    cgsetsprite(0)
    cgmakesprite(sprite.Resp_ItemPos,800,100, [1 0 0]) % ITEM POSITION
    cgtrncol(sprite.Resp_ItemPos,'r')
    cgsetsprite(sprite.Resp_ItemPos)
    cgrect(-300,0,150, 60, [0.3 0.3 0.3])
    cgrect(-95,0,150, 60, [0.3 0.3 0.3])
    cgrect(95,0,150, 60, [0.3 0.3 0.3])
    cgrect(300,0,150, 60, [0.3 0.3 0.3])
    cgtext('POSITION 1',-300,0)
    cgtext('POSITION 2',-95,0)
    cgtext('POSITION 3',95,0)
    cgtext('POSITION 4',300,0)
    cgsetsprite(0)
    
    % Item position options
    cgmakesprite(sprite.ItemPos,700,700, [1 0 0]) 
    cgtrncol(sprite.ItemPos,'r')
    cgsetsprite(sprite.ItemPos)
    cgrect(0,0,480, 300, w.back-0.5)
    cgpencol([0.8 0 0]);  cgfont('Helvetica',100)
    cgtext('1', -130, 75)
    cgtext('2', 130,75)
    cgtext('3', -130,-75)
    cgtext('4', 130,-75)
    cgpencol([0 0 0]);  cgfont('Helvetica',30)
    cgsetsprite(0)
    
    
end

%% EXECUTE MEMTEST

for trialnum=1:p.nTrials
    cgloadbmp(2,itemlist{data(trialnum, col.ItemStim),1},p.disp.ItemSize,p.disp.ItemSize)
    wk=[];
    
    % OLD/NEW?
    cgdrawsprite(2,0,25)
    cgdrawsprite(sprite.Resp_OldNew, 0,-250)
    cgtext('Is this item OLD or NEW?',0,200)
    cgtext('Have you seen this picture before?',0,-150)
    cgflip(w.back)
    clearkeys
    [wk.oldnewkey wk.time wk.n]= waitkeydown(inf,[p.keys.left p.keys.right]);
    switch wk.oldnewkey (1)
        case p.keys.left
            data(trialnum, col.OldNew)=1;
        case p.keys.right
            data(trialnum, col.OldNew)=2;
    end    
    
    % Remember/Know
    if data(trialnum, col.OldNew)==1;
        cgdrawsprite(2,0,25)
        cgdrawsprite(sprite.Resp_RemKnow, 0,-250)
        cgtext('Do you REMEMBER this item, ',0,250)
        cgtext('or do you just KNOW that it is old?',0,200)
        cgflip(w.back)
        clearkeys
        [wk.remknowkey wk.time wk.n]= waitkeydown(inf,[p.keys.left p.keys.right]);
        switch wk.remknowkey (1)
            case p.keys.left
                data(trialnum, col.RemKnow)=1;
            case p.keys.right
                data(trialnum, col.RemKnow)=2;
        end
    else
        data(trialnum, col.RemKnow)=999;
    end
    
    % SURE/GUESS
    cgdrawsprite(2,0,25)
    cgdrawsprite(sprite.Resp_SureGuess, 0,-250)
    if data(trialnum, col.OldNew)==1; cgtext('Are you SURE that this item is OLD?',0,225)
    else cgtext('Are you SURE that this item is NEW?',0,225)
    end
    cgflip(w.back)
    clearkeys
    [wk.sureguesskey wk.time wk.n]= waitkeydown(inf,[p.keys.left p.keys.right]);
    switch wk.sureguesskey (1)
        case p.keys.left
            data(trialnum, col.SureGuess)=1;
        case p.keys.right
            data(trialnum, col.SureGuess)=2;
    end
           
    % If Old, Item position?
    if data(trialnum, col.OldNew)==1;
        cgdrawsprite(sprite.ItemPos,0,0)
        cgdrawsprite(sprite.Resp_ItemPos, 0,-250)
        cgtext('Where on the screen did you see this item earlier?', 0,200)
        cgflip(w.back)
        clearkeys
        [wk.poskey wk.time wk.n]= waitkeydown(inf,[p.keys.k1 p.keys.k2 p.keys.k3 p.keys.k4 ]);
        switch wk.poskey(1)
            case p.keys.k1;  data(trialnum, col.Pos)=1;
            case p.keys.k2;  data(trialnum, col.Pos)=2;
            case p.keys.k3;  data(trialnum, col.Pos)=3;
            case p.keys.k4;  data(trialnum, col.Pos)=4;
        end
    end
    
    % Log responses as ROC
    switch [num2str(data(trialnum, col.OldNew)) num2str(data(trialnum, col.RemKnow)) num2str(data(trialnum, col.SureGuess))]
        case '111'; data(trialnum, col.Roc)=6;  % Old-Rem-Sure
        case '112';  data(trialnum, col.Roc)=5; % Old-Rem-Guess
        case '121'; data(trialnum, col.Roc)=4; % Old-Know-Sure
        case '122'; data(trialnum, col.Roc)=3; % Old-Know-Guess
        case '29992'; data(trialnum, col.Roc)=2; % New-Guess
        case '29991'; data(trialnum, col.Roc)=1; % New-Sure
        otherwise; data(trialnum, col.Roc)=999;
    end
    
    % Memory correct?
    if data(trialnum, col.OldNew)==data(trialnum, col.Item_OldNew);
        data(trialnum, col.CorrectRecognition)=1;
    else
        data(trialnum, col.CorrectRecognition)=0;
    end
    if data(trialnum, col.ItemPos)==data(trialnum, col.Pos);
        data(trialnum, col.CorrectPosition)=1;
    else
        data(trialnum, col.CorrectPosition)=0;
    end
    
    % Display
    disp(['ITEM #' num2str(data(trialnum, col.ItemStim)) '   (' itemlist{data(trialnum, col.ItemStim),1}(9:length(itemlist{data(trialnum, col.ItemStim),1})-4) ')'])
    disp(['Old/New: ' num2str(data(trialnum, col.OldNew))])
    disp(['Rem/Know: ' num2str(data(trialnum, col.RemKnow))])
    disp(['Sure/Guess: ' num2str(data(trialnum, col.SureGuess))])
    disp(['Roc: ' num2str(data(trialnum,col.Roc))]);
    disp(['Position: ' num2str(data(trialnum, col.Pos))])
    disp(['Correct Recognition? '  num2str(data(trialnum, col.CorrectRecognition))])
    disp(['Correct Position? '  num2str(data(trialnum, col.CorrectPosition))])
    disp('------------------------------------------------------------------------')
end

%%

cgflip(0.5,0.5,0.5)
cgtext('You have now completed',0,100)
cgtext('this part of the task',0,50)
cgtext('Please call the experimenter',0,-50)
t.end=cgflip(0.5,0.5,0.5);
p.log.end=clock;
waitkeydown(inf)
stop_cogent

%%

% Post-hoc calcs
outcomes.OldNewAccuracy=mean(data(:,col.CorrectRecognition));
outcomes.PositionAccuracy=mean(data(data(:,col.OldNew)==1,col.CorrectPosition));

% Save
mem.data=data;
mem.col=col;
mem.par=p;
mem.outcomes=outcomes;
save(['Data' filesep p.subject '_file_4Memtest.mat'], 'mem');
try % Transfer file to DeletedDaily 
    movetowhere='\\Asia\DeletedDaily\EL ExploreMem Data';
    if isdir(movetowhere)==0; mkdir(movetowhere);end
    copyfile(['Data' filesep p.subject '_file_4Memtest.mat'],  [movetowhere filesep p.subject '_file_4Memtest.mat']);
    w.warning=' ';
catch
    w.warning='WARNING: Data not transfered to Deleted Daily. Transfer files manually';
end

%
disp('--------------------------------------------------------------------------------------')
try % Duration 
    disp(['Duration: ' num2str((t.end-t.startcogent)/60) ' minutes'])
catch
end
disp(['Old/New Accuracy: ' num2str(outcomes.OldNewAccuracy)])
disp(['Position Accuracy: ' num2str(outcomes.PositionAccuracy)])
disp(w.warning)
disp('--------------------------------------------------------------------------------------')