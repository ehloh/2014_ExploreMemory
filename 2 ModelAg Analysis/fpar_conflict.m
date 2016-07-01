function [ data ] = fpar_conflict( data, col)
% [ data ] = fpar_conflict( data, col)
% Run through CONFLICT trials and mark RL parameters for each trial
%   
%   Input variables:                     EnvThreat (E), NTokens (N)
%                                               OutcomeMagnitude (optional)
%   Additional output variables:    pLoss (P), Entropy (U), VExplore (X), EntropyNTok (K)
%                                               (if OutcomeMagnitude available) OutcomeMean,  OutcomeVariance
%   
%   Columns for EnvThreat, NTokens, pLoss & Entropy must be specified. All
%           other variables only calculated if specified in variable col.
%
% --------------------------------------------------------------------------------------------------------------

% Execute to debug (get onsets): data=subjdata{s,2}.d_conflictcontrol(subjdata{s,2}.d_conflictcontrol(:,col.Task)==1,:);

% Setup
FixedLoss= - 8;
ExploreCost= - 2;
EntropyCorrection=0.00001;

%% Mark RL Variables (straightforward)

% (1) EnvThreat and NTokens: Convert to proper variable values if necessary
if sum(data(:, col.NTokens)==8)==0
    data(:,col.NTokens)=2*data(:,col.NTokens);
end
if sum(data(:,col.EnvThreat)<1)==0; 
    data(:,col.EnvThreat)=data(:,col.EnvThreat)/4;
end

% (2) pLoss
data(:,col.pLoss)=(data(:,col.EnvThreat)) .* (data(:,col.NTokens)/8);

% (3) Entropy
for i=1:size(data,1) % pLoss values  of 1 = corrected by -EntropyCorrection for entropy calculation
    
    % Convert if pLoss == 1 
    if data(i,col.pLoss)==1 
        data(i,col.pLoss)=data(i,col.pLoss)-EntropyCorrection;
    end
      
    % Entropy calculation:  -p(loggp)-(1-p)logg(1-p) 
    data(i,col.Entropy)= -(data(i,col.pLoss)).*log(data(i,col.pLoss)) -   (1-data(i,col.pLoss)).*log(1-data(i,col.pLoss));        
        
    % Convert pLoss back if it's been altered
    if data(i,col.pLoss)==1-EntropyCorrection;
        data(i,col.pLoss)=data(i,col.pLoss)+EntropyCorrection;
    end
    
end

% (4) Expected value
data(:,col.EV)=data(:,col.pLoss).*FixedLoss + (1-data(:,col.pLoss)).*data(:,col.NTokens);


%% Additional variables (that need thinking through)

% Value of exploration: Value gained from exploring
%       EV(AfterExplore) - EV(BeforeExplore)   [taking into account posterior probabilities after Exploring)
%       p(Accept/Reject after Explore) as binary, depending on which EV is bigger
if isfield(col, 'VExplore')
    
    % (a) Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
    %     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
    %                                      = (EnvThreat * NTokens) / (16 - EnvThreat*NTokens)
    ppActBomb=(data(:, col.EnvThreat).*data(:, col.NTokens)) ./  (16-(data(:, col.EnvThreat).*data(:, col.NTokens)));
    
    % (b) Likelihood of Accepting, if one Explores and sees no Bomb
    %   Cannot assume heuristic Accept/Reject deterministically from Explored-Info - will produce erroneously low V(Explore) values.
    %   For now, can assume deterministic post-Exploration behaviour:  deterministically  Accept/Reject
    %       depending on which has higher value. Can also calculate p(Accept) etc on the basis of the softmax betas.
    pAccGivNoBombSeen = (1-ppActBomb)   .* data(:, col.NTokens)>=   ppActBomb .* FixedLoss;
    
    % (c) EV(Explore)=   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
    %                       = p(See)*V(See)
    %                         +    p(~See) *  [       p(Accept|~See) * { posterior p(ActBomb) * Fixed Loss    +   (1- posterior p(ActBomb))   *  NTokens   }      +     p(Reject|~See) * 0      ]
    %                         +   Exploration Cost
    %     p(See Bomb) * V(See Bomb) = 0  ( Assume Reject, thus V(See)=0)
    %     p(Accept|~See): Deterministic/binary, if V(Accept)>V(Reject) based on posterior probabilities
    %     V(~See):   pAcceptGiven~See * [   posterior pActBomb * Fixed Loss + (1- posterior pActBomb )*NTokens  ]    + pRejectGiven~See * 0
    V_AccGivNoSee=    ppActBomb  .* FixedLoss           +            (1-ppActBomb ) .* data(:,col.NTokens);
    EV_NoSee= pAccGivNoBombSeen .*  V_AccGivNoSee  + (1-pAccGivNoBombSeen) .* 0;
    EV_Explore= 0 + (1-  data(:, col.pLoss).*0.5) .* EV_NoSee +  ExploreCost;
%     EV_Explore= 0 + data(:, col.pLoss).*0.5 .* EV_NoSee +  ExploreCost;  % This is wrong
    
    % (d) Value of Exploring =EV(Explore) - EV(NonExplore)
    %       EV(NonExplore)=EV(Accept) if >0, otherwise 0
    EV_NonExplore=data(:,col.EV).*(data(:, col.EV)>=0);
    data(:, col.VExplore)=EV_Explore-EV_NonExplore;
end

% EntropyNTok
if isfield(col, 'EntropyNTok')
    data(:, col.EntropyNTok)=data(:, col.Entropy).* data(:, col.NTokens);
end

% % Conflict: Uncertainty * (EVpositive - EVnegative)
% if isfield(col, 'Conflict')
%     data(:,col.Conflict)=data(:,col.Entropy).*((1-data(:,col.pLoss)).*data(:,col.NTokens) -   data(:,col.pLoss).*FixedLoss);
% end


% Theoretical Outcome variance/std dev and value according to mean-variance theory (D'acremont & Bossaerts 2008 CABN)
%     meanoutcome= Sum_forNoutcomes [ probability(outcome) * magnitude(outcome) ]
%     Unbiased sample variance= [1/(N-1)] .* Sum_forNoutcomes [    probability(outcome) * (outcome - meanoutcome)^2   ]
%     Std Deviation= sqrt( Unbiased sample variance )
%     * NOTE: Code here must agree with fcf_meanvar_quantities
OutcomeMean=data(:,col.EV);
StanDev=sqrt(   (1-data(:, col.pLoss)).*(data(:, col.NTokens)-OutcomeMean).^2 +  data(:, col.pLoss).*(FixedLoss -OutcomeMean).^2   );
if isfield(col,'StanDev')
    data(:, col.StanDev)=StanDev;
end
if isfield(col,'vMeanVar')
    data(:, col.vMeanVar)=  OutcomeMean - StanDev;   % V=mean - b*standev
end

% Binomial pActBomb variance
if isfield(col,'BinomVar')
    data(:, col.BinomVar)=data(:,col.pLoss).*(1-data(:,col.pLoss));
end

%% Outcome statistics (optional)

if isfield(col, 'OutcomeMagnitude')
    for e=1:4
        for n=2:2:8
            wo.dat=data(data(:, col.EnvThreat)==e/4 & data(:, col.NTokens)==n,col.OutcomeMagnitude);
            data(data(:, col.EnvThreat)==e & data(:, col.NTokens)==n, col.OutcomeMean)=mean(wo.dat);
            data(data(:, col.EnvThreat)==e & data(:, col.NTokens)==n, col.OutcomeVariance)=var(wo.dat);
            wo=[];
        end
    end
end

%% Plots to survey (request)

plot=0;
if plot==1
    % Settings for figure
    fontsize=15;
    fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms

    % Arrange data (all columns) by cell
    dmat=cell(4,4);
    for e=1:4
        for n=1:4
            dmat{5-e, n}=data(data(:,col.EnvThreat)== e/4 & data(:,col.NTokens)== n*2,:);
            
            % Individual variables
            et(5-e,n)=unique(dmat{5-e, n}(:, col.EnvThreat));
            nt(5-e,n)=unique(dmat{5-e, n}(:, col.NTokens));
            pl(5-e,n)=unique(dmat{5-e, n}(:, col.pLoss));
            en(5-e,n)=unique(dmat{5-e, n}(:, col.Entropy));
            ev(5-e,n)=unique(dmat{5-e, n}(:, col.EV));
            if isfield(col, 'VExplore');  vex(5-e,n)=unique(dmat{5-e, n}(1, col.VExplore)); end
            if isfield(col, 'EntropyNTok'); entok(5-e,n)=unique(dmat{5-e, n}(1, col.EntropyNTok)); end
            if isfield(col, 'EntropyNTok'); enev(5-e,n)=unique(dmat{5-e, n}(1, col.EntropyNTok)); end
            if isfield(col, 'Conflict'); cf(5-e,n)=unique(dmat{5-e, n}(1, col.Conflict)); end
            %
            if isfield(col, 'StanDev'); sd(5-e,n)=dmat{5-e, n}(1, col.StanDev); end
            if isfield(col, 'vMeanVar'); vmv(5-e,n)=dmat{5-e, n}(1, col.vMeanVar); end
            if isfield(col, 'BinomVar'); bnv(5-e,n)=dmat{5-e, n}(1, col.BinomVar); end
        end
    end
    close all hidden; figure('Name', 'Conflict task RL variables', 'NumberTitle', 'off', 'Position', [100 100 1400 800], 'Color', 'w');
    
    % EnvThreat
    subplot(2,5,1); imagesc(et,[0 1]); title('Env Threat')
    set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar; 
    
    % NTokens
    subplot(2,5,2);  imagesc(nt,[2 12]);  title('N Tokens')
    set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar; 
    
    % pLoss
    subplot(2,5,3); imagesc(pl ,[0 1]); title('pLoss')
    set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar; 
    
    % Entropy
    subplot(2,5,4); imagesc(en); title('Entropy')
    set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar; 
    
    % Conflict
%     subplot(2,5,5); imagesc(cf); title('Conflict [Entropy*(EVpos-EVneg)]')
%     set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar; 
   
    % #################################################################
     
    % EV
    subplot(2,5,5); imagesc(ev,[min(min(ev)) max(max(ev))]); title('EV')
    set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar; 

    % VExplore
    if isfield(col, 'VExplore');  
        subplot(2,5,6); imagesc(vex); title('VExplore (pprob)')
        set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar;
        % subplot(2,5,7); imagesc(vex, [min(min(vex)) max(max(vex))]); title('VExplore (min=0)')
        %         set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar;
    end
    
    % EntropyNTok
    if isfield(col, 'EntropyNTok');  
        subplot(2,5,7); imagesc(entok,[min(min(entok)) max(max(entok))]); title('Entropy * NTok')
        set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar;
    end
    
    % #################################################################
     
   
    % #################################################################
     
    % Binomial variance
    if isfield(col, 'BinomVar');  
         subplot(2,5,8); imagesc(bnv); title('Binomial variance')
         set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar;
    end
    
   % Theoretical standev + mean-variance value
    if isfield(col, 'StanDev');  
         subplot(2,5,9); imagesc(sd); title('Standard Dev')
         set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar;
    end
    if isfield(col, 'vMeanVar');  
         subplot(2,5,10); imagesc(vmv); title('Mean-Var value')
         set(gca,'YTick', 1:4, 'YTickLabel', fliplr({'1/4' '2/4' '3/4' '1'}),'XTick', 1:4, 'XTickLabel', 2:2:8, 'FontSize', fontsize); axis square; colorbar;
    end
    
    
    
end  % for plot

end


