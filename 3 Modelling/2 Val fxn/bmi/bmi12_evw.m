function [ V ] = bmi12_evw(x, modelinput)
% Inoptimal value function (post-prob multiplier) + distortion of V(Accept|~See) 
%      + Fixed exploration bonus (e) + Variable exploration bonus (w) x VExplore (v)
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
%                                        x(2)= Unoptimal posterior probability distortion   (m)
%                                        x(3)= V(Accept|~See) distortion (i)
%                                        x(4)= Fixed exploration bonus (e)
%                                        x(5)= Variable exploration bonus (w)
%-------------------------------------------------------------------------------------------------------------------------------

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

ExploreBonusFixed=x(:,4);
ExploreBonusVar=x(:,5);

%% EV calculation using basic value function

% Get parameters in model
modname=mfilename; modname(strfind(modname, '_')-2:strfind(modname, '_'))=[];

[ V ] = bmi01(x, {modname data, fPar,col});

% Add variable Exploration bonus x VExplore
nTrials=size(data,1); nPP=size(x,1); 
V(:,:,3)=V(:,:,3)  +   repmat(ExploreBonusFixed, [1 nTrials]);
V(:,:,3)=V(:,:,3)  +   repmat(ExploreBonusVar, [1 nTrials]).*repmat(data(:,col.VExplore)', [nPP 1]);

end

