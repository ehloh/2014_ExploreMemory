% Run Env-Learning session (trial parameters set up previous script)
%   To alter trial parameters, do so in 1st script
clear all; close all hidden; clc

for o1=1:1 % Documentation
    % This script sets up columns 1-6; executed in other scripts during tasks themselves
    %         Col 1:       TrialType
    %         Col 2:       EnvThreat (1=Low, 2=High)
    %         Col 3:       NToken pairs (1 or 3)
    %         Col 4:       Bomb present?
    %         Col 5:       Activated bomb present?
    %         Col 6:       Trial no.
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
    col.Resp1=8;
    col.RT1=9;
    col.Resp2=10;
    col.RT2=11;
    col.BombExplored=12;
    col.Trialvalid=13;
    col.OutcomeMag=14;
    col.ItemStim=15;
    col.ItemType=16;
end
for o1=1:1 % Execution
    testing=1;
    
    if testing==1
        p.subject=input('Subject ID:   ', 's');
        w.fullscreen=1;
    else
        p.subject='t1';
        w.fullscreen=0;
    end
end
for o1=1:1 % General settings + load subject's trial parameters
    pp=load(['Data' filesep p.subject '_file_1taskpar.mat']);
    p=pp.Learn.par;
    rdata=pp.Learn.learnpar;
    col=pp.Learn.col;
    p.nTrials=size(rdata,1);
    rdata(:,col.Trialnum)=rand(p.nTrials,1);
    rdata=sortrows(rdata, col.Trialnum); rdata(:,col.Trialnum)=1:p.nTrials;
    %
    rdata(:,col.Resp1)=1; % LearnEnv: forced accept
    
    % Misc
    w.line=50;
    w.back=p.disp.settokencolor;
    w.res=2;
    rand('state',sum(100*clock));
    
end

% 
p.keys.monitor=71;

% COGENT
config_display(w.fullscreen,w.res, [0 0 0], [1 1 1], 'Helvetica', 40, 9,0, 0); % marker 5g start
config_keyboard;
start_cogent % if using testing laptops, w.res in config dis must be 3! also, w.lineoptions =-330;
% cgtext('The computer will now give you', 0, w.line*4)
% cgtext('instructions for this task. ', 0, w.line*3)
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
    sprite.ResponseOptions=5;
    sprite.TokenActive=6;
    sprite.Explored_Bomb=7;
    sprite.Explored_NotBomb=8;
    sprite.Outcome_WinGreen=9;
    sprite.Outcome_LoseRed=10;
    sprite.EnvThreat_Offers=[11 12 13 14]; % Empty-token offers for all levels of EnvThreat
    
    % (A) Make sprites for each EnvThreat type (all empty tokens)
    for i=1:p.design.EnvThreat_Levels
        cgmakesprite(sprite.EnvThreat_Offers(i),1000,700, p.disp.contingencycolours{i,1})
        cgsetsprite(sprite.EnvThreat_Offers(i))
        for j=1:p.design.NSpaces_TokenPairs
            cgpencol(p.disp.settokencolor)
            cgpenwid(8)
            cgellipse(p.disp.tokenpos(j),p.disp.tokenY_row1,p.disp.tokensize,p.disp.tokensize)
            cgellipse(p.disp.tokenpos(j),p.disp.tokenY_row2,p.disp.tokensize,p.disp.tokensize)
        end
        cgsetsprite(0)
    end
    
    % (B) Create all other template sprites
    for o2=1:1
        
        % [Offer] Offered token
        cgmakesprite(sprite.TokenActive,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
        cgtrncol(sprite.TokenActive,'n')
        cgsetsprite(sprite.TokenActive)
        cgpencol([p.disp.settokencolor])
        cgellipse(0,0,p.disp.tokensize,p.disp.tokensize,'f')
        
        % [Explore] Bomb token
        cgmakesprite(sprite.Explored_Bomb,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
        cgtrncol(sprite.Explored_Bomb,'n')
        cgsetsprite(sprite.Explored_Bomb)
        cgpenwid(10); cgpencol(p.disp.color_bomb)
        cgellipse(0,0,p.disp.tokensize,p.disp.tokensize)
        
        % [Explore] Not a bomb
        cgmakesprite(sprite.Explored_NotBomb,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
        cgtrncol(sprite.Explored_NotBomb,'n')
        cgsetsprite(sprite.Explored_NotBomb)
        cgpenwid(10); cgpencol(p.disp.color_nobomb)
        cgellipse(0,0,p.disp.tokensize,p.disp.tokensize)
        
        % [Outcome] Win token
        cgmakesprite(sprite.Outcome_WinGreen,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
        cgtrncol(sprite.Outcome_WinGreen,'n')
        cgsetsprite(sprite.Outcome_WinGreen)
        cgpencol([p.disp.color_nobomb])
        cgellipse(0,0,p.disp.tokensize+8,p.disp.tokensize+8,'f')
        
        % [Outcome] Loss token
        cgmakesprite(sprite.Outcome_LoseRed,p.disp.tokensize+10,p.disp.tokensize+10,[0 0 0])
        cgtrncol(sprite.Outcome_LoseRed,'n')
        cgsetsprite(sprite.Outcome_LoseRed)
        cgpencol([p.disp.color_bomb])
        cgellipse(0,0,p.disp.tokensize+8,p.disp.tokensize+8,'f')
        
        % Response options
        cgmakesprite(sprite.ResponseOptions,700,100, [1 0 0])
        cgtrncol(5,'r')
        cgsetsprite(5)
        cgrect(-230,0,170, 60, [0.3 0.3 0.3])
        cgrect(0,0,170, 60, [0.3 0.3 0.3])
        cgrect(230,0,170, 60, [0.3 0.3 0.3])
        cgpencol([0 0 0])
        cgfont('Helvetica',30)
        cgtext('ACCEPT',-230,0)
        cgtext('REJECT',0,0)
        cgtext('EXPLORE',230,0)
        cgfont('Helvetica',40)
        cgsetsprite(0)
    end
    
    % (C) Test display (optional)
    w.testdisplay=0;
    if w.testdisplay==1
        
        % Test colours
        w.width=100;
        for i=1:6; % p.design.EnvThreat_Levels
            w.col=p.disp.contingencycolours{i};
            cgmakesprite(20+i,w.width,600,w.col)
            cgdrawsprite(20+i,-(w.width*3)-40+i*w.width,0)
        end
        cgflip(p.disp.settokencolor)
        
        % Test offer displays
        for i=1:p.design.EnvThreat_Levels
            cgdrawsprite(sprite.EnvThreat_Offers(i),0,0)
            cgflip(0,0.1,0)
            waitkeydown(inf)
        end
        
        % Test tokens
        cgpencol(p.disp.settokencolor)
        w.sprites2disp=fieldnames(sprite);
        cgtext('START TOKENS',0,0); cgflip(0,0.1,0)
        waitkeydown(inf)
        for i=1:length(w.sprites2disp)
            eval(['w.nsprite=sprite.' w.sprites2disp{i} ';'])
            cgdrawsprite(sprite.EnvThreat_Offers(1),0,0)
            cgdrawsprite(w.nsprite, p.disp.tokenpos(randi(length(p.disp.tokenpos),1)),p.disp.tokenY_row1)
            cgdrawsprite(w.nsprite,p.disp.tokenpos(randi(length(p.disp.tokenpos),1)),p.disp.tokenY_row2)
            cgflip(0,0,0)
            waitkeydown(inf)
        end
    end
    cgpencol(0,0,0);
end

%% EXECUTE TASK

w.block=1;
cgtext('+',0,0)
countmonitor=0;

% error

for trialnum=1: p.nTrials
    
    % Verify parameters-display match-up
    disp('--------------------------')
    disp(['Trial #: ' num2str(trialnum)])
    disp(['Offer: Colour ' num2str(rdata(trialnum,col.EnvThreat)) ', ' num2str(2*rdata(trialnum,col.NTokenPairs)) ' activated'])
    disp(['Active Bomb=' num2str(rdata(trialnum,col.BombActivated))])
    %     wtt.goon=waitkeydown(inf);
    
    for o1=1:1 % Set up displays & details for this trial
        wt=[]; wk=[]; % Resests
        
        % Shuffle positions of display
        wt.pos=[p.disp.tokenpos' zeros(length(p.disp.tokenpos ),1)];
        wt.startpos=randi(p.design.NSpaces_TokenPairs+1-rdata(trialnum,col.NTokenPairs),1);
        for i=wt.startpos:(wt.startpos+rdata(trialnum,col.NTokenPairs)-1)
            wt.pos(i,2)=1;
        end
        wt.posX(1:rdata(trialnum,col.NTokenPairs),1)=wt.pos((wt.pos(:,2)==1),1); % Positions in which to mark activated tokens
        
    end
    rdata(trialnum,col.Trialnum)=trialnum;
    wt.fixate=cgflip(w.back-0.5);
 
    % [FIRST OFFER ONLY] -----------
    cgdrawsprite(sprite.EnvThreat_Offers(rdata(trialnum,col.EnvThreat)),0,0)
    for j=1:rdata(trialnum,col.NTokenPairs)
        cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row1);
        cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row2);
    end
        cgtext('+',0,0)
%     cgdrawsprite(sprite.ResponseOptions,0,-250)
    %
    waituntil(wt.fixate*1000+500+rand*500)
    clearkeys
    wt.offer1only=cgflip(w.back);
    wt.offer1item=wt.offer1only;
    
    % [FIRST OFFER + Item] -----------
%     cgdrawsprite(sprite.EnvThreat_Offers(rdata(trialnum,col.EnvThreat)),0,0)
%     for j=1:rdata(trialnum,col.NTokenPairs)
%         cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row1);
%         cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row2);
%     end
%     cgdrawsprite(sprite.ResponseOptions,0,-250)
%     %
%     waituntil(wt.offer1only*1000+p.disp.tOffer1Only)
%     clearkeys
%     wt.offer1item=cgflip(w.back);
    rdata(trialnum,col.Resp1)=1; % Assumed Acceptance
    rdata(trialnum,col.Trialvalid)=1;
    
    % [EXPLORE/2ND OFFER: WILL NEVER HAPPEN] -------------
    if rdata(trialnum,col.Resp1)==3
        cgdrawsprite(sprite.EnvThreat_Offers(rdata(trialnum,col.EnvThreat)),0,0)
        for j=1:rdata(trialnum,col.NTokenPairs)
            cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row2);
        end
        cgdrawsprite(sprite.ResponseOptions,0,-250)
        %
        cgrect(230,-250,120, 30, [0.1 0.1 0.1])
        eval(['wt.revealed_Y=p.disp.tokenY_row' num2str(randi(2,1)) ';'])  % Reveal status of half
        for j=1:rdata(trialnum,col.NTokenPairs) % Mark Explored-NotBombs
            cgdrawsprite(sprite.Explored_NotBomb,wt.posX(j),wt.revealed_Y);
        end
        rdata(trialnum, col.BombExplored)=randi(2)-1;
        %
        if rdata(trialnum, col.BombActivated)==1 && rdata(trialnum, col.BombExplored)==1 % Bomb activated + Bomb explored
            wt.bombposx=wt.posX(randi(size(wt.posX,1)));
            rdata(trialnum,col.BombExplored)=1;
            cgdrawsprite(sprite.Explored_Bomb,  wt.bombposx,  wt.revealed_Y);
        else
            rdata(trialnum,col.BombExplored)=0;
        end
        cgtext(['Information cost: - ' num2str(p.design.ExplorationCost*10) ' p'],0,-130)
        cgtext('+',0,0)
        waituntil(wt.offer1item*1000+p.disp.tOffer1Item)
        clearkeys
        wt.offer2=cgflip(w.back);
        [wk.key2press wk.key2time wk.key2n]=waitkeydown(p.disp.tOffer2, [wt.key1 wt.key2]);
        if wk.key2n==0 % No response
            cgdrawsprite(sprite.EnvThreat_Offers(rdata(trialnum,col.EnvThreat)),0,0)
            for j=1:rdata(trialnum,col.NTokenPairs)
                cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row1);
                cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row2);
            end
            cgdrawsprite(sprite.ResponseOptions,0,-250)
            cgrect(230,-200,120, 30, [0.1 0.1 0.1])
            for j=1:rdata(trialnum,col.NTokenPairs) % Mark new information
                cgdrawsprite(sprite.Explored_NotBomb,wt.posX(j),wt.revealed_Y);
            end
            if rdata(trialnum, col.BombActivated)==1 && rdata(trialnum, col.BombExplored)==1 % Bomb activated + Bomb explored
                cgdrawsprite(sprite.Explored_Bomb,  wt.bombposx,  wt.revealed_Y);
            end
            cgtext(['Information cost: - ' num2str(p.task.explorationcost*10) ' p'],0,-130)
            %
            cgtext('NO RESPONSE',0,160)
            waituntil(wt.offer2*1000+p.disp.tOffer2)
            wt.outcome=cgflip(w.back);
            rdata(trialnum,col.Trialvalid)=0;
        else
            rdata(trialnum,col.Trialvalid)=1;
            switch wk.key2press(1) % log responses
                case wt.key1
                    rdata(trialnum,col.Resp2)=1;
                case wt.key2
                    rdata(trialnum,col.Resp2)=2;
            end
            rdata(trialnum,col.RT2)= wk.key2time(1)-wt.offer2*1000;
        end
    else
        wt.offer2=wt.offer1item;
        rdata(trialnum,col.BombExplored)=999;
    end
    
    % [Calculate outcome]-------
    if rdata(trialnum,col.Trialvalid)==1
        if rdata(trialnum,col.Resp1)==3 % Explore?
            wt.explorecost=p.design.ExplorationCost;
            wt.infocost=['Information cost: - ' num2str(p.design.ExplorationCost*10) ' p'];
            wt.respcol=col.Resp2;
        else
            wt.explorecost=0;
            wt.infocost=' ';
            wt.respcol=col.Resp1;
        end
        if rdata(trialnum,wt.respcol)==2 % Reject: nothing
            rdata(trialnum,col.OutcomeMag)=0-wt.explorecost;
        elseif rdata(trialnum,wt.respcol)==1 && rdata(trialnum,col.BombActivated)==1 % Accept +Bomb
            rdata(trialnum,col.OutcomeMag)=p.design.FixedLoss-wt.explorecost;
        elseif rdata(trialnum,wt.respcol)==1 && rdata(trialnum,col.BombActivated)==0 % Accept +No Bomb
            rdata(trialnum,col.OutcomeMag)=2*rdata(trialnum,col.NTokenPairs)-wt.explorecost;
        end
        if rdata(trialnum,col.OutcomeMag)>0 % Display colours
            wt.won=['+ ' num2str((rdata(trialnum,col.OutcomeMag))*10) ' p'];
            wt.tokenspritenum=sprite.Outcome_WinGreen;
        elseif rdata(trialnum,col.OutcomeMag)==0
            if rdata(trialnum,wt.respcol)==1
                wt.won='+ 0 p';
                wt.tokenspritenum=sprite.TokenActive;
            else
                wt.won='?';
                wt.tokenspritenum=sprite.TokenActive;
            end
        else
            wt.won=[num2str(rdata(trialnum,col.OutcomeMag)*10) ' p'];
            wt.tokenspritenum=sprite.Outcome_LoseRed;
            if rdata(trialnum,wt.respcol)==2
                wt.won='?';
                wt.tokenspritenum=sprite.TokenActive;
            end
        end
    end
    
    % Display outcome
    if rdata(trialnum,col.Trialvalid)==1
        cgdrawsprite(sprite.EnvThreat_Offers(rdata(trialnum,col.EnvThreat)),0,0)
        for j=1:rdata(trialnum,col.NTokenPairs)
            cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row2);
        end
        %             cgdrawsprite(sprite.ResponseOptions,0,-250)
        %
        for j=1:rdata(trialnum,col.NTokenPairs) % win/loss tokens
            cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row2);
        end
        cgfont('Helvetica',70)
        cgtext(wt.won,0,160)
        cgfont('Helvetica',40)
        cgtext(wt.infocost,0,-130)
        cgtext('+',0,0)
        waituntil(wt.offer2*1000+p.disp.tOffer2)
        wt.outcome=cgflip(w.back);
    end
    
    % Monitoring task?
    if rdata(trialnum,col.MonitorTask)==1
        cgdrawsprite(sprite.EnvThreat_Offers(rdata(trialnum,col.EnvThreat)),0,0)
        for j=1:rdata(trialnum,col.NTokenPairs)
            cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(sprite.TokenActive,wt.posX(j),p.disp.tokenY_row2);
        end
        for j=1:rdata(trialnum,col.NTokenPairs) % win/loss tokens
            cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row1);
            cgdrawsprite(wt.tokenspritenum,wt.posX(j),p.disp.tokenY_row2);
        end
        cgfont('Helvetica',70)
        cgtext(wt.won,0,160)
        cgfont('Helvetica',40)
        cgtext(wt.infocost,0,-130)
        %        
        cgpencol(1,1,1);
        cgtext('+',0,0)
        cgpencol(0,0,0);
        %
        waituntil(wt.outcome*1000+700+randi(300))
        wt.showmonitor=cgflip(w.back);
        clearkeys
        [wk.monitorkeypress wk.monitorkeytime wk.monitorkeyn]=waitkeydown(500, p.keys.monitor); % Half a second to respond
        countmonitor=countmonitor+wk.monitorkeyn;
        disp(['Monitor task caught?   ' num2str(wk.monitorkeyn)])
        waituntil(wt.showmonitor*1000+1000)
    end
    
    
    % Breaks
    if trialnum==floor(size(rdata,1)*w.block/4)
        wait(2000)
        cgflip(0.5,0.5,0.5)
        cgtext('Please take a break if you like',0,100)
        cgtext(['You have completed ' num2str(w.block) '/4 of this task'],0,0)
        cgtext('Press any key to continue',0,-100)
        w.block=w.block+1;
        t.break=cgflip(0.5,0.5,0.5);
        waitkeydown(inf)
    end
    
    cgtext('+',0,0)
    waituntil(wt.outcome*1000+p.disp.tOutcome)
    
    % Email researcher if almost done
    if trialnum==size(rdata,1)-15
        try
            w.time=clock;
            f_sendemail('learnreward',['[ExploreSimple] (' num2str(w.time(4)) ':', num2str(w.time(5)) 'hrs) ' w.subjname ' has 1.5 min left (' w.taskstage ')'], ['Session almost complete: Conflict task'])
        end
    end
    
    % Verify choices recorded
    disp(['Response valid: ' num2str(rdata(trialnum,col.Trialvalid))])
    disp(['Response 1:' num2str(rdata(trialnum,col.Resp1)) '  (1=Accept, 2=Reject, 3=Explore)'])
    disp(['Response 2:' num2str(rdata(trialnum, col.Resp2)) '  (1=Accept, 2=Reject, 3=Explore)'])
    disp(['Outcome: ' num2str(rdata(trialnum,col.OutcomeMag))])
    %         wtt.goon=waitkeydown(inf,60);
    
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

% Rate EnvThreat probabilities (experimenter to execute)
disp('##################### END ######################')
disp('###############################################')
disp('###############################################'); disp(' '); disp(' ');
disp('################### RATING TASK ###################');disp(' '); 
disp('EXPERIMENTER: Enter subject''s estimated p(Bomb) for each colour, as decimals or fractions.'); 
disp('                         (rating probability of bomb - NOT probability of an ACTIVATED bomb)'); disp(' ');disp(' ');

rate_envthreat=cell(p.design.EnvThreat_Levels, 3); % [1] Colour [2] EnvThreatLevel (1=Low) [3] Rated probability
rate_envthreat(:,[1 2])=[p.disp.contingencycolours(1: p.design.EnvThreat_Levels,2) num2cell(1: p.design.EnvThreat_Levels)'];
for i=1: p.design.EnvThreat_Levels
    disp(['Q' num2str(i) ': How likely is there to be a bomb (activated or inactivated) on ' rate_envthreat{i,1} ' trials?'])
    rate_envthreat{i,3}=input(['      (pBomb) on ' rate_envthreat{i,1} ' trials:   ']);
end
w.rate=sortrows(rate_envthreat,3);
outcomes.EnvThreat_Accuracy=mean(cell2mat(w.rate(:,2))==cell2mat(rate_envthreat(:,2)));
disp('###############################################')

%%

% Post-hoc calcs
outcomes.monitorok=countmonitor/sum(rdata(:,col.MonitorTask));
outcomes.n_omissions=sum(rdata(:,col.Trialvalid)==0);
outcomes.n_explorations=sum(rdata(:,col.Resp1)==3);
outcomes.net=sum(rdata((rdata(:,col.Trialvalid)==1),col.OutcomeMag));

% Save
p.log.date=date;
cF.rdata=rdata;
cF.rdata=rdata(rdata(:,col.Trialvalid)==1,:);
cF.par=p;
cF.col=col;
cF.outcomes=outcomes;
cF.rate_envthreat=rate_envthreat;
[cF.behaviour] = f_plotindividualbehaviour(p.design,rdata, col); % Plot behaviour
save(['Data' filesep p.subject '_file_2LearnEnv.mat'], 'cF');
try % Transfer file to DeletedDaily 
    movetowhere='\\Asia\DeletedDaily\EL ExploreSimplified Data';
    if isdir(movetowhere)==0; mkdir(movetowhere);end
    copyfile(['Data' filesep p.subject '_file_2LearnEnv.mat'],  [movetowhere filesep p.subject '_file_2LearnEnv.mat']);
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
disp(['Score: ' num2str(outcomes.net) ' (' num2str(outcomes.net/  sum(rdata(rdata(:,col.BombActivated)==0, col.NTokenPairs))*2)   '  % of theoretical max)'])
disp(['# Omissions: ' num2str(outcomes.n_omissions)])
disp(['# Explored: ' num2str(outcomes.n_explorations)])
disp(['EnvThreat acc: ' num2str(outcomes.EnvThreat_Accuracy)   '    (p.disp.contingencycolours)']);
disp(['Monitoring task accuracy:  ' num2str(outcomes.monitorok*100)  ' %']);
disp(w.warning)
disp('--------------------------------------------------------------------------------------')