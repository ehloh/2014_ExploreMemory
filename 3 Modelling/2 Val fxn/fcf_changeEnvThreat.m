function [ outvars] = fcf_changeEnvThreat(EnvThreat, NTok, miscvar)
% [ outvars] = fcf_changeEnvThreat(EnvThreat, NTok, miscvar)
% Given EnvThreat + NTok vectors, output all other parameter vectors
%   
%   Input variables:                     EnvThreat (E), NTokens (N)
%                                               OutcomeMagnitude (optional)
%   Additional output variables:    pLoss (P), Entropy (U), VExplore (X), EntropyNTok (K)
%                                               (if OutcomeMagnitude available) OutcomeMean,  OutcomeVariance
%   
%   Columns for EnvThreat, NTokens, pLoss & Entropy must be specified. All
%           other variables only calculated if specified in variable col.

try  FixedLoss=miscvar.FixedLoss; catch  FixedLoss= -8; end
try  ExploreCost=miscvar.ExploreCost; catch  ExploreCost= -2; end
EntropyCorrection=0.00001;

%%

% (1) EnvThreat and NTokens: Convert to proper variable values if necessary
if sum(NTok==8)==0; NTok=2*NTok;
end
if sum(EnvThreat<1)==0;  EnvThreat=EnvThreat/4;
end

% (2) pLoss
pLoss=(EnvThreat).*(NTok./8);

% (3) Entropy
pLossEn=pLoss;  pLossEn(pLossEn==1)=  pLossEn(pLossEn==1)- EntropyCorrection;
Entropy= -(pLossEn).*log(pLossEn) -   (1-pLossEn).*log(1-pLossEn);

% (4) EntropyNTok
EntropyNTok=Entropy.*NTok;

% (5) EV
EV=pLoss.*FixedLoss + (1-pLoss).*NTok;

% (6) VExplore   ##############################
[ ppActBomb] = f_posteriorprob( EnvThreat, NTok);

% (b) Likelihood of Accepting, if one Explores and sees no Bomb
%   Cannot assume heuristic Accept/Reject deterministically from Explored-Info - will produce erroneously low V(Explore) values.
%   For now, can assume deterministic post-Exploration behaviour:  deterministically  Accept/Reject
%       depending on which has higher value. Can also calculate p(Accept) etc on the basis of the softmax betas.
pAccGivNoBombSeen = (1-ppActBomb)   .* NTok>=   ppActBomb .* FixedLoss;

% (c) EV(Explore)=   p(See)*V(See)     +     p(~See)* V(~See)    +   Exploration Cost
%                       = p(See)*V(See)
%                         +    p(~See) *  [       p(Accept|~See) * { posterior p(ActBomb) * Fixed Loss    +   (1- posterior p(ActBomb))   *  NTokens   }      +     p(Reject|~See) * 0      ]
%                         +   Exploration Cost
%     p(See Bomb) * V(See Bomb) = 0  ( Assume Reject, thus V(See)=0)
%     p(Accept|~See): Deterministic/binary, if V(Accept)>V(Reject) based on posterior probabilities
%     V(~See):   pAcceptGiven~See * [   posterior pActBomb * Fixed Loss + (1- posterior pActBomb )*NTokens  ]    + pRejectGiven~See * 0
V_AccGivNoSee=    ppActBomb  .* FixedLoss           +            (1-ppActBomb ) .* NTok;
EV_NoSee= pAccGivNoBombSeen .*  V_AccGivNoSee  + (1-pAccGivNoBombSeen) .* 0;
EV_Explore= 0 + (1-pLoss.*0.5) .* EV_NoSee +  ExploreCost;

% (d) Value of Exploring =EV(Explore) - EV(NonExplore)
%       EV(NonExplore)=EV(Accept) if >0, otherwise 0
EV_NonExplore= EV.*(EV>=0);
VExplore=EV_Explore-EV_NonExplore;

%% Output

outvars.EnvThreat =EnvThreat ;
outvars.NTok =NTok ;
outvars.NTokens =NTok ;
outvars.pLoss =pLoss;
outvars.Entropy=Entropy;
outvars.EntropyNTok =EntropyNTok;
outvars.VExplore=VExplore;
outvars.EV=EV;

end

