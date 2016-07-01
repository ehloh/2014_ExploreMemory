function [ V ] = bma01(x, modelinput)
% Model's value function: [ V ] = bma01(x, modelinput)
%   Value function with inoptimal posterior probability calculation (m), distortion of V(Accept|~See)
%
%  - Inputs
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[]    data   fPar   col}
%                                   data:   Data from all trials     (row=trialnum, col=fixed params in col)
%                                   col:      Column specifications for data
%                                   fPar:    Fixed parameters
%
%  -  Output 
%        V:        Value associated with Accept (1), Reject (2), Explore (3)
%                        (row=parampoint, col=trialnum, 3rd dimension=Accept/Reject/Explore)
%
% ---------------------------------------------------------------------------------------------------------------------------------
% Free parameters 'x':        x(1)= softmaxbeta        (Applied outside of value fxn)
%                                        x(2)= Unoptimal posterior probability distortion   (m)
%                                        x(3)= V(Accept|~See) distortion (a)
% ---------------------------------------------------------------------------------------------------------------------------------

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};
beta=20./(1+exp(-x(:,1)));          % 0 < beta <20
UnOptim= 20./(1+exp(x(:,2))) ;       % 0 < m < 20
Distort_VAccGivNoSee=x(:,3);

nPP=size(x,1);  % No. of parameter points (unique combinations of parameter values for all free parameters covered by grid)

% % % Simplifications for debugging.
% data=sortrows(data,[col.EnvThreat col.NTokens]);  data=data(1:18:648,:); 
% cnames=cellstr([repmat('e', 36,1) num2str(data(:,col.EnvThreat)*6) repmat('n', 36,1) num2str(data(:,col.NTokens)/2)])';

%% Set up details

% Data
data_cF=data(data(:, col.Task)==1,:); nT_cF = size(data_cF,1);
order_cF=data(data(:, col.Task)==1,col.Trialnum); 
data_ct=data(data(:, col.Task)==2,:);  nT_ct= size(data_ct,1);
order_ct=data(data(:, col.Task)==2,col.Trialnum); 

% Fixed parameters that sometimes act as free parameters) - convert to array if necessary
if length(fPar.cF_FL)==1; fPar.cF_FL=repmat(fPar.cF_FL,[nPP 1]); end
if length(fPar.ct_FL)==1; fPar.ct_FL=repmat(fPar.ct_FL,[nPP 1]); end
if length(fPar.cF_EC)==1; fPar.cF_EC=repmat(fPar.cF_EC,[nPP 1]); end
if length(fPar.ct_EC)==1; fPar.ct_EC=repmat(fPar.ct_EC,[nPP 1]); end


%% Conflict task - Basic EV function

V_cF=nan(nPP, nT_cF ,3);

% [1] V(Accept) & V(Reject)
V_cF(:,:,1)= repmat(data_cF(:,col.pLoss)',[nPP 1]) .* repmat(fPar.cF_FL, [1 nT_cF]) + repmat(  (   (1- data_cF(:,col.pLoss)).*data_cF(:,col.NTokens))', [nPP 1]);
V_cF(:,:,2)=0;

% [2a] Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
%     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
%                                      = (EnvThreat * NTokens) / (24 - EnvThreat*NTokens)
ppActBomb=(data_cF(:, col.EnvThreat).*data_cF(:, col.NTokens)) ./  (24-(data_cF(:, col.EnvThreat).*data_cF(:, col.NTokens)));
ppActBomb=repmat(ppActBomb', [nPP 1]);
ppABscaled=power(ppActBomb, repmat(UnOptim, [1 nT_cF]));

% [2b] Likelihood of Accepting, if one Explores and sees no Bomb
%   p(Accept) is calculated using the softmax (with existing betas) and the posterior values 
%       (posVChoice, V Accept/Reject after Exploring)
%   Cannot assume that subjects will heuristically Accept if Exploring finds no Activated Bombs. 
%       Calculating EV(Explore) must weight the value by the probability of Accepting/Rejecting.
%       Otherwise, V(Exploring) will sometimes given values that are lower than the EC, which 
%       would not occur because it is more optimal to Reject such offers. Previously, had assumed deterministic  
%       post-Exploration behaviour:  deterministically  Accept/Reject depending on which has higher value. 
% pAccGivNoBombSeen = (1-ppActBomb)   .* repmat(data_cF(:, col.NTokens)', [nPP 1])  +  ppActBomb.* repmat(fPar.cF_FL, [1 nT_cF])>=0;
posVChoice=nan(size(V_cF));
posVChoice(:,:,1)=ppABscaled.* repmat(fPar.cF_FL, [1 nT_cF]) + (1-ppABscaled).* repmat(data_cF(:,col.NTokens)', [nPP 1]);
posVChoice=posVChoice(:); posVChoice(posVChoice==0)=eps; posVChoice=reshape(posVChoice, size(V_cF)); % remove V(Accept)=0.
posVChoice(:,:,2)=0;
pAccGivNoBombSeen=exp((repmat(beta, [1 nT_cF]).*posVChoice(:,:,1)))   ./   (   exp(repmat(beta, [1 nT_cF]).*posVChoice(:,:,1))  +   exp(repmat(beta, [1 nT_cF]).*posVChoice(:,:,2))    );

% [3] V(Explore)=   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
%                   = p(See)*V(See)     
%                         +    p(~See) *  [       p(Accept|~See) * { posterior p(ActBomb) * Fixed Loss    +   (1- posterior p(ActBomb))   *  NTokens   }      +     p(Reject|~See) * 0      ]
%                         +   Exploration Cost
%     p(See Bomb) * V(See Bomb) = 0  ( Assume Reject, thus V(See)=0)
%     p(Accept|~See): Deterministic/binary, if V(Accept)>V(Reject) based on posterior probabilities
%     V(~See):   pAcceptGiven~See * [   posterior pActBomb * Fixed Loss + (1- posterior pActBomb )*NTokens  ]    + pRejectGiven~See * 0
V_AccGivNoSee=   ppABscaled  .* repmat(fPar.cF_FL,[1 nT_cF]) +            (1-ppABscaled) .* repmat(data_cF(:,col.NTokens)', [nPP 1]);
V_AccGivNoSee=V_AccGivNoSee + Distort_VAccGivNoSee;
EV_NoSee= pAccGivNoBombSeen .*  V_AccGivNoSee  + (1-pAccGivNoBombSeen) .* 0;
V_cF(:,:,3)= 0 +  repmat((data_cF(:, col.pLoss).*0.5)', [nPP 1])  .* EV_NoSee +  repmat(fPar.cF_EC, [1 nT_cF]);

%% Control task - Basic EV function

V_ct=nan(nPP, nT_ct ,3);

% [1] V(NoBomb) & V(Bomb)
V_ct(:,:,1)= repmat(((1- data_ct(:,col.pLoss)).*data_ct(:,col.NTokens))', [nPP 1]);
V_ct(:,:,2)= repmat((data_ct(:,col.pLoss).*data_ct(:,col.NTokens))', [nPP 1]);

% [2a] Posterior probability of there being an Activated Bomb, if Exploring reveals no bomb (i.e. ~See)
%     posterior p(ActBomb)= p(~See | UnseenActBomb) * p(UnseenActBomb) / (1 -EnvThreat*N/24)
%                                      = (EnvThreat * NTokens) / (24 - EnvThreat*NTokens)
%   Kept in vector format, like data: row=trialnum
ppActBomb=(data_ct(:, col.EnvThreat).*data_ct(:, col.NTokens))  ./  (24-(data_ct(:, col.EnvThreat).*data_ct(:, col.NTokens)));
ppActBomb=repmat(ppActBomb', [nPP 1]);
ppABscaled=power(ppActBomb,repmat(UnOptim, [1 nT_ct]));

% [2b] Likelihood of Accepting, if one Explores and sees no Bomb
%   Kept in vector format, like data: row=trialnum
posVChoice=nan(size(V_ct));
posVChoice(:,:,1)= repmat(data_ct(:, col.NTokens)', [nPP 1])  .* (1-ppABscaled);
posVChoice(:,:,2)= repmat(data_ct(:, col.NTokens)', [nPP 1])  .* ppABscaled;
pAccGivNoBombSeen=exp((repmat(beta, [1 nT_ct]).*posVChoice(:,:,1)))   ./   (   exp(repmat(beta, [1 nT_ct]).*posVChoice(:,:,1))  +   exp(repmat(beta, [1 nT_ct]).*posVChoice(:,:,2))    );

% [3] V(Explore) =   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
VAccGivNoSee=  pAccGivNoBombSeen.*posVChoice(:,:,1);
VAccGivNoSee= VAccGivNoSee + Distort_VAccGivNoSee;
VRejGivNoSee=(1-pAccGivNoBombSeen).*posVChoice(:,:,2);
EVGivenNotSee =VAccGivNoSee+VRejGivNoSee;
% EVGivenNotSee =pAccGivNoBombSeen.*posVChoice(:,:,1)+(1-pAccGivNoBombSeen).*posVChoice(:,:,2);
V_ct(:,:,3)=  repmat((data_ct(:, col.pLoss)*0.5)', [nPP 1])  .* repmat(data_ct(:, col.NTokens)', [nPP 1])+  repmat((1-(data_ct(:, col.pLoss)*0.5))', [nPP 1]).*EVGivenNotSee  + repmat(fPar.ct_EC, [1 nT_ct]);

%% Combine values for both tasks

Vv=[V_cF V_ct];
% V=Vv; 
% 
% Preserve trial order (as experienced by subject - other functions expect this)
Vv(nPP+1,:,:)=repmat([order_cF; order_ct]',[1 1 3]);
Vv(:,:,1)=sortrows(Vv(:,:,1)',nPP+1)';
Vv(:,:,2)=sortrows(Vv(:,:,2)',nPP+1)';
Vv(:,:,3)=sortrows(Vv(:,:,3)',nPP+1)';

V=nan(nPP,nT_cF+nT_ct,3);
V(:,:,1)=Vv(1:nPP,:,1);
V(:,:,2)=Vv(1:nPP,:,2);
V(:,:,3)=Vv(1:nPP,:,3);


end
