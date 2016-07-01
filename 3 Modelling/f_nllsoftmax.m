function [nll,pch]=f_nllsoftmax(x, modelinput, varargin)
% [nll,pch]=f_nllsoftmax(x, modelinput)
% Calculate neg log-likelihood associated with model
%
%   Inputs
%       x                   Array of parameter points (all unique combinations of parameter values for 
%                               each free parameter) (row= param-point, col=parameter)
%
%       modelinput    Value function for model, data, fixed param & column specifications (for data)
%                               i.e. {ModelValueFxn      data   fixedpars     col}
%
%       Z                   Priors for hierarchical fit (Z.mu= mean, Z.nui=variance)
%
%       doprior           Apply penalty for hierarchical fit?
% ------------------------------------------------------------------------------------- 

% Additional non-standard inputs?
if isempty(varargin)==0;
    Z=varargin{1};doprior=varargin{2};
else doprior=0;
end

% beta=exp(x(:,1)); % Beta is constrained here - all other parameters constrained in value function
beta=20./(1+exp(-x(:,1)));
ValueFxn=modelinput{1}; data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

%%

nTrials=size(data,1);
nParPoints=size(x,1);

% Apply model's value function to derive values 
%       (V: row=param-point, col=trialnum, 3rd dimension=Accep/Reject/Explore)
eval(['V=' ValueFxn '(x, {[] data, fPar,col});']) 

% If necessary, correct values to allow for non-nan exponential values
checkb=exp(   repmat(beta, [1 nTrials 3]) .* V );
if sum(isinf(checkb(:)))~=0
    Vcorrection=max(V(:));
    V=V-Vcorrection;
else Vcorrection=0;
end

% Calculate p(Choice) using softmax rule
%           p(Choice)=exp(beta*V(Choice))/(exp(beta*V(Option 1))+exp(beta*V(Option 2))+exp(beta*V(Option 3)));
b2 = repmat(beta, [1 nTrials]);
base=(exp(b2.*squeeze(V(:,:,1)))  +   exp(b2.*squeeze(V(:,:,2)))    +   exp(b2.*squeeze(V(:,:,3)))); % point in parameter space * trial * ac/re/ex
if sum(base==0)~=0; base(base==0)=eps;   disp(['   Correcting for softmax base ==0 (beta=' num2str(beta) ')']); end % Correct for values close to 0 (occurs with extreme betas)
actualChoice = data(:,col.Choice);
Vchosen = repmat((actualChoice==1)', [nParPoints 1]).*V(:,:,1)+repmat((actualChoice==2)', [nParPoints 1]).*V(:,:,2)+repmat((actualChoice==3)', [nParPoints 1]).*V(:,:,3);
pch=exp(b2.*Vchosen)./base;

% Restore original value magnitudes 
V=V+Vcorrection;

% Calculate models negative log-likelihood
nll=  - sum(log(pch),2);

%% [Hierarchical fit] Penalize extreme parameter values 
%   Group parameter distribution (assumed gaussian) is used to constrain individual parameters (by applying a nLL penalty to extreme values)
%   Make sure Z.mu is in the same space as x

if doprior
    % Calculating the penalty, as a function of mean (mu) and variance (nui) of parameter distributions
    LL_penalty=    -1/2 * (x-Z.mu)   *  Z.nui  * (x-Z.mu)' - 1/2*log(2*pi/det(Z.nui));
    nll = nll  - LL_penalty;

    %   Quentin's original: lp = -1/2 * (x-Z.mu)'*Z.nui*(x-Z.mu) - 1/2*log(2*pi/det(Z.nui));
    %                                l  = -l  - sum(lp);
    
%     % Penalty calculation #2
%     LL_penalty = log(mvnpdf(x, Z.mu, Z.nui));
%     nll  = nll  - LL_penalty;
%     if mvnpdf(x, Z.mu, Z.nui)==0; 
%         disp('nLL penalty is infinite. Probability density at current parval =0')
%         for pn=1:length(Z.mu)
%             subplot(length(Z.mu),1,pn)
%             a=normrnd(  Z.mu(pn), Z.nu(pn), 10000,1 );
%             hist(a); hold on;  scatter(x(pn), 500,'r')
%             title(['mu=' num2str(Z.mu(pn)) ', nu=' num2str( Z.nu(pn))])
%         end
%         nll=  - sum(log(pch),2)  -  log(eps);  % Redo penalty calculation
%     end
end

   

end


