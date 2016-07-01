function [ V ] = bp10_fow(x, modelinput)
% Value function + Subjective value of Loss + Variable exploration bonus (w) x EntropyNTok  (o)
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
%                                        x(2)= softmax epsilon  (p)
%                                        x(3)= Subjective loss value (f)
%                                        x(4)= Variable exploration bonus (w)
%-------------------------------------------------------------------------------------------------------------------------------

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

fPar.cF_FL=x(:,3);
ExploreBonus=x(:,4);

%% EV calculation using basic value function

% Get parameters in model
modname=mfilename; modname(strfind(modname, '_')-2:strfind(modname, '_'))=[];

[ V ] = bp01(x, {modname data, fPar,col});

% Add variable Exploration bonus x EntropyNTok (Entropy x NTokens)
nTrials=size(data,1); nPP=size(x,1); 
%
[ cfvars] = fcf_meanvar_quantities(data(data(:, col.Task)==1, col.pLoss), data(data(:, col.Task)==1, col.NTokens), fPar);
[ ctvars] = fct_meanvar_quantities(data(data(:, col.Task)==2, col.pLoss), data(data(:, col.Task)==2, col.NTokens), fPar);
stddev=nan(size(data,1),1); stddev(data(:,col.Task) ==1)=cfvars.stddev; stddev(data(:,col.Task) ==2)=ctvars.stddev;
V(:,:,3)=V(:,:,3)  +  repmat(ExploreBonus, [1 nTrials]).*repmat(  stddev', [nPP 1]);

end

