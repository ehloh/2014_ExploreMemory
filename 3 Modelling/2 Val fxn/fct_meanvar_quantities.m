function [ ov] = fct_meanvar_quantities(prob, ntok, invar)
% [ outvars] = f_meanvar_quantities(task, prob, ntok, invar)
% Calculate Std Dev and Value (according to Mean-Variance theory; see D'acremont & Bossaerts 2008 CABN)
% 
%  [Inputs]:
%         task, prob, ntok:        Vectors (same length)
%         invar.b                       b, uncertainty-sensitivity (see value output)
%                                                 (if absent, b=1)
%             
% [Outputs]:
%         stddev:                         Standard deviation (in value calc)
%         meanoutcome:             EV, mean outcome (in value calc)
%         value:                          V= meanoutcome - b*StdDev
%
% NOTE: Code here must agree with fpar_control
%
% #############################################

try  mv_b=invar.b; catch  mv_b=1; end
EntropyCorrection=0.00001;

%%

% Entropy calculation:  -p*(logp)  -   (1-p)*[ log(1-p) ]
entropy= -prob.*(log(prob))  - (1-prob).*(log(1-prob));
entropy(find(prob==1))= -(1-EntropyCorrection).*(log(1-EntropyCorrection))  - (1-(1-EntropyCorrection)).*(log(1-(1-EntropyCorrection)));

% Mean outcome/EV = p(Correct)*V(Correct) + p(Wrong)*V(Wrong)
mo=(1-entropy).*ntok  + entropy.*0;
stddev=sqrt(   (1-entropy).*(ntok-mo).^2 +  entropy.*(0 -mo).^2   );
v=mo - mv_b.*stddev;

%% Output

ov.stddev=stddev;
ov.meanoutcome=mo;
ov.value=v;

end

