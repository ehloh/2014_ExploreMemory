function [ V ] = b13_eow(x, modelinput)
%   Basic value function + Fixed exploration bonus (e) + Variable exploration bonus (w) x EntropyNTok (o)
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
%                                        x(2)= Fixed exploration bonus (e)
%                                        x(3)= Variable exploration bonus (w)
%-------------------------------------------------------------------------------------------------------------------------------

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

ExploreBonusFixed=x(:,2);
ExploreBonusVar=x(:,3);

%% EV calculation using basic value function

% Get parameters in model
modname=mfilename; modname(strfind(modname, '_')-2:strfind(modname, '_'))=[];

[ V ] = b01(x, {modname data, fPar,col});

% Add variable Exploration bonus x Variance [pLoss*(1-pLoss)]
nTrials=size(data,1); nPP=size(x,1); 
V(:,:,3)=V(:,:,3)  +   repmat(ExploreBonusFixed, [1 nTrials]);
V(:,:,3)=V(:,:,3)  +  repmat(ExploreBonusVar, [1 nTrials]).*repmat( (    data(:, col.pLoss).*(1-data(:, col.pLoss)) )', [nPP 1]);

end