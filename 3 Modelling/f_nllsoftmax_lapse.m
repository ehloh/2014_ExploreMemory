function [nll,pch]=f_nllsoftmax_lapse(x, modelinput, varargin)
% [nll,pch]=f_nllsoftmax_lapse(x, modelinput)
% Calculate neg log-likelihood associated with model
%
%   Inputs
%       x                   Array of parameter points (all unique combinations of parameter values for 
%                               each free parameter) (row= param-point, col=parameter)
%                               [1st param=beta, 2nd param=epsilon]
%
%       modelinput      Value function for model, data, fixed param & column specifications (for data)
%                               i.e. {ModelValueFxn      data   fixedpars     col}
%
%       Z                   Priors for hierarchical fit (Z.mu= mean, Z.nui=variance)
%
%       doprior           Apply penalty for hierarchical fit?
% ------------------------------------------------------------------------------------- 

% Additional non-standard inputs?
if isempty(varargin)==0;
    Z=varargin{1}; doprior=varargin{2};
else doprior=0;
end

% Softmax parameters constrained here - all others constrained in value function
beta=20./(1+exp(-x(:,1)));  % 0 < beta <20
epsilon=1./(3+3.*exp(-x(:,2)));   %  0 < epsilon < 1/3 (sigmoid function)

ValueFxn=modelinput{1}; data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

%%

nTrials=size(data,1);
nParPoints=size(x,1);

% Apply model's value function to derive values 
%       (V: row=param-point, col=trialnum, 3rd dimension=Accep/Reject/Explore
eval(['V=' ValueFxn '(x, {[] data, fPar,col});']) 

% If necessary, correct values to allow for non-nan exponential values
checkb=exp(   repmat(beta, [1 nTrials 3]) .* V );
if sum(isinf(checkb(:)))~=0
    Vcorrection=max(V(:));
    V=V-Vcorrection;
else Vcorrection=0;
end


% %
% badparvals=f_transpar({'b'; 'p'; 'm'; 'a'; 'f'; 'e'; 'w'}, x, 'to');
% ts=[data(:, [col.EnvThreat col.NTokens  col.pLoss])  shiftdim(V,1) ]; vv= shiftdim(V,1);
% sum((ts(:,3).*badparvals(5)+ (1- ts(:,3)) .* ts(:,2)) - ts(:,4))  % All V(Accept) ok?
% p=ts(:,3); pp= (p./2)./(1-p./2);
% vANS1= power(pp, badparvals(3)) .*  badparvals(5)  +  (1-power(pp, badparvals(3))).*ts(:, 2);
% vANS=1./(1+badparvals(4).*(-vANS1));
% vExp= p.*0.5.*0  +    (1-p.*0.5).*(vANS)  -2  + badparvals(6)+ badparvals(7)*data(:, col.EntropyNTok)  ;
% discrep=vExp-ts(:,6);


% Calculate p(Choice) associated with Value functions, using softmax rule + epsilon
%       p(Choice)=exp(beta*V(Choice))/(exp(beta*V(Option 1))+exp(beta*V(Option 2))+exp(beta*V(Option 3)));
b2 = repmat(beta, [1 nTrials]);
base=(exp(b2.*squeeze(V(:,:,1)))  +   exp(b2.*squeeze(V(:,:,2)))    +   exp(b2.*squeeze(V(:,:,3)))); % point in parameter space * trial * ac/re/ex
if sum(base(:)==0)~=0;   bb=base(:); bb(bb==0)=eps; base=reshape(bb, size(b2)); disp(['   Correcting for softmax base ==0 (beta=' num2str(beta) ')']); end  % Correct for values close to 0 (occurs with extreme betas)
actualChoice = data(:,col.Choice);
Vchosen = repmat((actualChoice==1)', [nParPoints 1]).*V(:,:,1)+repmat((actualChoice==2)', [nParPoints 1]).*V(:,:,2)+repmat((actualChoice==3)', [nParPoints 1]).*V(:,:,3);
pch=exp(b2.*Vchosen)./base;

% Restore original value magnitudes 
V=V+Vcorrection;

% Episilon: Model accepts some stochasticity in choice
%       p(Choice)=  epsilon + (1-3*epsilon)* p(Choice)      [p(Choice) according to existing softmax calculation]
%                       =  epsilon + (1-3*epsilon).*( exp(beta*V(Choice))/(exp(beta*V(Option 1))+exp(beta*V(Option 2))+exp(beta*V(Option 3)))   )
%       The logic of it is that choice is random some amount of the time (reflected in constant epsilon added to pChoice)
%       Epsilon is thus bound by chance p(Choice) - in this case, 1/3
pch= repmat(epsilon, [1 nTrials]) + (1-3.*repmat(epsilon, [1 nTrials])).*pch;

% Calculate models negative log-likelihood
nll=  - sum(log(pch),2);


%% [Hierarchical fit] Penalize extreme parameter values 
%   Group parameter distribution (assumed gaussian) is used to constrain individual parameters (by applying a nLL penalty to extreme values)

if doprior
    
    % Calculating the penalty, as a function of mean (mu) and variance (nui) of parameter distributions
    LL_penalty=    -1/2 * (x-Z.mu)   *  Z.nui  * (x-Z.mu)' - 1/2*log(2*pi/det(Z.nui));
    nll  = nll  - LL_penalty;
%     
%     LL_penalty=    -1/2 * (x-Z.mu)   *  Z.nui  * (x-Z.mu)' - 1/2*log(2*pi/det(Z.nui));
%     nll  = nll  - LL_penalty;
    
    
    %   Quentin's original: lp = -1/2 * (x-Z.mu)'*Z.nui*(x-Z.mu) - 1/2*log(2*pi/det(Z.nui));
    %                                l  = -l  - sum(lp);
    
%     
%      LL_penalty = log(mvnpdf(x, Z.mu, Z.nui));
%      nll  = nll  - LL_penalty;
end

end


