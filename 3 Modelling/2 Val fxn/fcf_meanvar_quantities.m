function [ ov] = fcf_meanvar_quantities(prob, ntok, invar)
% [ outvars] = f_meanvar_quantities(task, prob, ntok, invar)
% Calculate Std Dev and Value (according to Mean-Variance theory; see D'acremont & Bossaerts 2008 CABN)
% 
%  [Inputs]:
%         task, prob, ntok:        Vectors (same length)
%         invar.b                       b, uncertainty-sensitivity (see value output)
%                                                 (if absent, b=1)
%         invar.FixedLoss           Subjective value of loss
%                                                 (if absent, FixedLoss= -12)
%             
% [Outputs]:
%         stddev:                         Standard deviation (in value calc)
%         meanoutcome:             EV, mean outcome (in value calc)
%         value:                          V= meanoutcome - b*StdDev
%
% NOTE: Code here must agree with fpar_conflict
%
% #############################################

try  FixedLoss=invar.FixedLoss; catch  FixedLoss= -8; end
try  mv_b=invar.b; catch  mv_b=1; end

%%

mo= prob.*FixedLoss + (1-prob).*ntok;
stddev=sqrt(   (1-prob).*(ntok-mo).^2 +  prob.*(FixedLoss -mo).^2   );
v=mo - mv_b.*stddev;

%% Output

ov.stddev=stddev;
ov.meanoutcome=mo;
ov.value=v;

end

