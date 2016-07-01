function [ V ] = bi06_ow(x, modelinputs)
%  Distortion of V(Accept|~See)  + Variable exploration bonus (w) x EntropyNTok (o)
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
%                                        x(2)= V(Accept|~See) distortion (i)
%                                        x(3)= Variable exploration bonus (w)
%-------------------------------------------------------------------------------------------------------------------------------

data=modelinputs{2}; fPar=modelinputs{3}; col=modelinputs{4};

ExploreBonus=x(:,3);

%% EV calculation using basic value function

% Get parameters in model
modname=mfilename; modname(strfind(modname, '_')-2:strfind(modname, '_'))=[];

[ V ] = bi01(x, {modname data, fPar,col});

% Add variable Exploration bonus x Variance [pLoss*(1-pLoss)]
nTrials=size(data,1); nPP=size(x,1); 
V(:,:,3)=V(:,:,3)  +  repmat(ExploreBonus, [1 nTrials]).*repmat( (    data(:, col.pLoss).*(1-data(:, col.pLoss)) )', [nPP 1]);

end

