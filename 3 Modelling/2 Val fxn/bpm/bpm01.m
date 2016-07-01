function [ V ] = bpm01(x, modelinput)
%
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[]    data   fPar   col}
%        V:                   Value associated with Accept (1), Reject (2), Explore (3)
%                                   (row parampoint, col=trialnum, 3rd dimension=Accept/Reject/Explore)
%
% ---------------------------------------------------------------------------------------------------------------------------------

% Get model parameters if absent
if isempty(modelinput{1})==1;
    modname=mfilename; modname(strfind(modname, '0'):end)=[];
    modelinput{1}=modname;
end

% Remove p parameter - others act like they would in the corresponding model without p
%           epsilon is NOT fed into basic normative value function
x=x(:, [1 strfind(modelinput{1}, 'p')+1:end]);  
modelinput{1}(strfind(modelinput{1}, 'p'))=[];

[ V ] = bm01(x, modelinput);

end



