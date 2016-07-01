function [ V ] = bm15_fevw(x, modelinput)
% Inoptimal value function (post-prob multiplier) + Subjective value of Loss + Fixed exploration bonus (e) 
%                                   + Variable exploration bonus (w) x VExplore (v)
%
%  - Inputs
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[] data   fPar   col}
%                                   data:   Data from all trials     (row=trialnum, col=fixed params in col)
%                                   col:      Column specifications for data
%                                   fPar:    Fixed parameters
%
%  -  Output 
%        V:        Value associated with Accept (1), Reject (2), Explore (3)
%                        (row=parampoint, col=trialnum, 3rd dimension=Accept/Reject/Explore)
%
% ---------------------------------------------------------------------------------------------------------------------------------
% Free parameters 'x':        x(1)= softmax beta        (Applied outside of value fxn)
%                                        x(2)= Unoptimal posterior probability multipler  (m)
%                                        x(3)= Subjective loss value (f)
%                                        x(4)= Fixed exploration bonus (e)
%                                        x(5)= Variable exploration bonus (w)
%-------------------------------------------------------------------------------------------------------------------------------

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

EnvThreatDistort= 20./(1+exp(-x(:,3))) ;       % 0 < j < 20
fPar.cF_FL=x(:,4);
ExploreBonusFixed=x(:,5);
ExploreBonusVar=x(:,6);

%% EV calculation using basic value function

% Get parameters in model
modname=mfilename; modname(strfind(modname, '_')-2:strfind(modname, '_'))=[];

[ V ] = bpj01(x, {modname data, fPar,col});

% Update subjective pars for variable exploration bonus
mv.FixedLoss=fPar.cF_FL;
data(:, col.EnvThreat)=power(data(:, col.EnvThreat),EnvThreatDistort);
dcf=data(data(:,col.Task)==1, [col.EnvThreat col.NTokens]);  % cf
[ov] = fcf_changeEnvThreat(dcf(:,1), dcf(:,2), mv);
dct=data(data(:,col.Task)==2, [col.EnvThreat col.NTokens]);  % ct
[ovv] = fct_changeEnvThreat(dct(:,1), dct(:,2), mv);
data(data(:,col.Task)==1, [col.Entropy col.EntropyNTok col.VExplore]) = [ov.Entropy ov.EntropyNTok ov.VExplore];  % Change only variable explore bonus pars
data(data(:,col.Task)==2, [col.Entropy col.EntropyNTok col.VExplore]) = [ovv.Entropy ovv.EntropyNTok ovv.VExplore];  

% Add variable exploration bonus
nTrials=size(data,1); nPP=size(x,1); 
V(:,:,3)=V(:,:,3)  +   repmat(ExploreBonusFixed, [1 nTrials]);
V(:,:,3)=V(:,:,3)  +   repmat(ExploreBonusVar, [1 nTrials]).*repmat(data(:,col.VExplore)', [nPP 1]);

end

