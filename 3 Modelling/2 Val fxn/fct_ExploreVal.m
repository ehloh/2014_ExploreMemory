function [ data] = fcf_ExploreVal(data, col, subvar)
% [ outvars] = fcf_ExploreVal(EnvThreat, NTok, miscvar)
% Calculate/output miscellaneous quantities associated with exploration
% Code here should match the value functions in the models themselves 
%   
%   Input variables:                     
%             subvar   [structure-stored param names and values, e.g. subvar.f=-8]
%             data       [mandatory cols: EnvThreat, NTokens]
%
%   Output variables (request in col. *)
%             VExplore                Value of exploration, before exploration
%             EV_expbomb         V(Gamble) after one explores and sees a bomb
%             EV_nobomb          V(Gamble) after one explores and sees NO bomb
%             
% ####################################################

% Load and check fitted param values 
if isfield(subvar, 'ExploreCost')==1; ExploreCost=subvar.ExploreCost; else ExploreCost=-2; end
if isfield(subvar, 'vw')==1; error('VExplore (as an exploration bonus multiplier has not been set up yet, for calculating sub-quantities in VExplore'); end
if isfield(subvar, 'j')==1; error('Has j param been implented properly yet? CHECK. '); end
if isfield(subvar, 'm')==1; m=subvar.m; else m=1; end
if isfield(subvar, 'i')==1; i=subvar.i; else i=0; end
% if isfield(subvar, 'f')==1; FixedLoss=subvar.f; else FixedLoss= -8; end
if isfield(subvar, 'e')==1; e=subvar.e; else e=0; end
    
EntropyCorrection=0.00001;

%%

if isfield(subvar, 'j')==1;
    error('Haven''t worked out j param yet!')
    % Distort EnvThreat + everything knocks on
%     EnvThreat=power(data_cF(:,col.EnvThreat), EnvThreatDistort);
%     [ov] = fcf_changeEnvThreat(EnvThreat, data_cF(:,col.NTokens), miscvar);    
end

% [1] Assemble + check existing values (envthreat, ntokens, ploss, entropy)
if sum(data(:, col.NTokens)==8)==0; data(:, col.NTokens)=2*data(:, col.NTokens); end
if sum(data(:, col.EnvThreat)<1)==0;  data(:, col.EnvThreat)=data(:, col.EnvThreat)/4; end
if isfield(col, 'pLoss')==0; data(:, col.pLoss)=(data(:, col.EnvThreat)).*(data(:, col.NTok)./8); end
if isfield(col, 'Entropy')==0; %  Entropy: -p*(logp)-(1-p)*(log(1-p))
    data(:, col.Entropy)=-data(:, col.pLoss).*log(data(:, col.pLoss))  -(1-data(:, col.pLoss)).*log((1-data(:, col.pLoss)));
    data(data(:, col.pLoss)==1, col.Entropy)= -(1-EntropyCorrection)*log(1-EntropyCorrection) -  (1-(1-EntropyCorrection))* log(1-(1-EntropyCorrection));
end
if isfield(col, 'EntropyNTok')==0; data(:, col.EntropyNTok)=data(:, col.Entropy).*data(:, col.NTokens); end


if isfield(subvar, 'e') & isfield(subvar, 'i')
    i=i-e;
end


%% Calculate VExplore-related quantities

% [2a] Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
[ ppActBomb] = f_posteriorprob( EnvThreat, NTok);
ppABscaled=power(ppActBomb, m);

% [2b] Likelihood of Accepting, if one Explores and sees no Bomb
posVacc= (1-ppABscaled) .* data(:, col.NTokens);
posVacc=posVacc+i;
%
posVrej= ppABscaled .* data(:, col.NTokens);
pAccGivNoBombSeen=posVacc>posVrej;
EVGivenNotSee = pAccGivNoBombSeen.* posVacc   +   (1-pAccGivNoBombSeen).* posVrej;  

% [3] V(Explore) =   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
VExplore =  (  data(:, col.pLoss).*0.5) .* data(:, col.NTokens)  +   (1- data(:, col.pLoss).*0.5) .* EVGivenNotSee   + ExploreCost;

% [4] Add exploration bonuses
VExplore =VExplore  +e; 
if isfield(subvar,'w')==1
    if isfield(subvar,'ow')==1;                  VExplore = VExplore  +  subvar.w.*data(:, col.EntropyNTok);
    elseif isfield(subvar,'uw')==1;            VExplore = VExplore  +  subvar.w.*data(:, col.Entropy);
    elseif isfield(subvar,'vw')==1;             error('Ugh set up vw param in VExplore-related quants')
    end
end


%% Output requested variables

if isfield(col, 'vSee')
    data(:, col.vSee)=data(:, col.NTokens); 
end

if isfield(col, 'vNoSee')
    data(:, col.vNoSee)=EVGivenNotSee;
end


end


