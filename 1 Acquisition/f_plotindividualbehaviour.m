function [b] = f_plotindividualbehaviour(p,data, col)
% Generate probability heat-plots for individual subject's behaviour
% [b] = f_plotindividualbehaviour(p,data, subjname)
%
% 
% Specify:
%
%         p.EnvThreat_Levels
%         p.NTokenPairs_Levels
%         p.EnvThreat_pBomb
%                     - if empty, assume pBomb = increments of 1/p.EnvThreat_Levels,
%                         with maximal EnvThreat p(Bomb)=1. - i.e. (1:p.EnvThreat_Levels)/p.EnvThreat_Levels;
%         p.NTokenPairs_Npairs
%                     - if empty, assume 1:1: n_levels
%         col.EnvThreat
%         col.NTokenPairs
%         col.OutcomeMag
%         col.Resp1
%         col.RT1
%
% Output - plots & variables:
% 
%         b.accept
%         b.reject
%         b.explore
%         b.RT
%         b.outcomemean
%         b.outcomevar
%
% Note: An earlier version of this exists with different (less flexible) specifications (25/9/13)
% -------------------------------------------------------------------------------------------------

% Execute to debug: data=rdata; p=p.design;

for o1=1:1 % Checks + Fill in the blanks for unspecified parameters
    if isfield(p,'EnvThreat_pBomb')==0
        p.EnvThreat_pBomb=(1:p.EnvThreat_Levels)/p.EnvThreat_Levels;
    end
    if isfield(p, 'NTokenPairs_Npairs')==0
        p.NTokenPairs_Npairs=1:1:p.NTokenPairs_Levels;
    end
end

%% Generate data matrices

b.accept=nan*ones(p.EnvThreat_Levels,p.NTokenPairs_Levels);  % Resets
b.reject=nan*ones(p.EnvThreat_Levels,p.NTokenPairs_Levels);
b.explore=nan*ones(p.EnvThreat_Levels,p.NTokenPairs_Levels);
b.RT=nan*ones(p.EnvThreat_Levels,p.NTokenPairs_Levels);
b.outcomemean=nan*ones(p.EnvThreat_Levels,p.NTokenPairs_Levels);
b.outcomevar=nan*ones(p.EnvThreat_Levels,p.NTokenPairs_Levels);

for e=1:p.EnvThreat_Levels 
    for n=1:p.NTokenPairs_Levels 
        ws=[];
        ws.choice=data(data(:,col.EnvThreat)==e & data(:,col.NTokenPairs)==p.NTokenPairs_Npairs(n) , col.Resp1);
        ws.RT=data(data(:,col.EnvThreat)==e & data(:,col.NTokenPairs)==p.NTokenPairs_Npairs(n) , col.RT1);
        ws.outcome=data(data(:,col.EnvThreat)==e & data(:,col.NTokenPairs)==p.NTokenPairs_Npairs(n) , col.OutcomeMag);
        %
        b.accept(p.EnvThreat_Levels+1- e,  n)=sum(ws.choice(:)==1)/size(ws.choice,1);
        b.reject(p.EnvThreat_Levels+1- e,  n)=sum(ws.choice(:)==2)/size(ws.choice,1);
        b.explore(p.EnvThreat_Levels+1- e,  n)=sum(ws.choice(:)==3)/size(ws.choice,1);
        b.RT(p.EnvThreat_Levels+1- e,  n)=mean(ws.RT);
        b.outcomemean(p.EnvThreat_Levels+1-e,  n)=mean(ws.outcome);
        b.outcomevar(p.EnvThreat_Levels+1- e,  n)=var(ws.outcome);
    end
end

%% Plot

figure('Name', 'Behavioural performance' ,'NumberTitle','off','Position',[715,220,900,600]);
set(gcf,'Color',[1 1 1])
for e=1:p.EnvThreat_Levels
    ylabels{p.EnvThreat_Levels+1-e}=strtrim(rats(p.EnvThreat_pBomb(e)));
end
for n=1:p.NTokenPairs_Levels
    xlabels{n}=num2str(2*p.NTokenPairs_Npairs(n));
end
    
    subplot(2,3,1) % Accept
    imagesc(b.accept, [0 1])
    axis('square'); title('% Accept'); colorbar
    ylabel('EnvThreat');
    set(gca,'YTick',1:p.EnvThreat_Levels)
    set(gca, 'YTickLabel', ylabels)    
    xlabel('NTokens');
    set(gca,'XTick',1:p.NTokenPairs_Levels)
    set(gca, 'XTickLabel',xlabels)
    
    subplot(2,3,2) % Reject
    imagesc(b.reject, [0 1])
    axis('square'); title('% Reject'); colorbar
    ylabel('EnvThreat');
    set(gca,'YTick',1:p.EnvThreat_Levels)
    set(gca, 'YTickLabel', ylabels)    
    xlabel('NTokens');
    set(gca,'XTick',1:p.NTokenPairs_Levels)
    set(gca, 'XTickLabel',xlabels)
    
    subplot(2,3,3) % Explore
    imagesc(b.explore, [0 1])
    axis('square'); title('% Explore'); colorbar
    ylabel('EnvThreat');
    set(gca,'YTick',1:p.EnvThreat_Levels)
    set(gca, 'YTickLabel', ylabels)    
    xlabel('NTokens');
    set(gca,'XTick',1:p.NTokenPairs_Levels)
    set(gca, 'XTickLabel',xlabels)
    
    subplot(2,3,4) % Mean RT
    imagesc(b.RT)
    axis('square'); title('Mean RT'); colorbar
    ylabel('EnvThreat');
    set(gca,'YTick',1:p.EnvThreat_Levels)
    set(gca, 'YTickLabel', ylabels)    
    xlabel('NTokens');
    set(gca,'XTick',1:p.NTokenPairs_Levels)
    set(gca, 'XTickLabel',xlabels)
    
    subplot(2,3,5) % Outcome mean
    imagesc(b.outcomemean)
    axis('square'); title('Outcome Mean'); colorbar
    ylabel('EnvThreat');
    set(gca,'YTick',1:p.EnvThreat_Levels)
    set(gca, 'YTickLabel', ylabels)    
    xlabel('NTokens');
    set(gca,'XTick',1:p.NTokenPairs_Levels)
    set(gca, 'XTickLabel',xlabels)
    
    subplot(2,3,5) % Outcome variance
    imagesc(b.outcomevar)
    axis('square'); title('Outcome Variance'); colorbar
    ylabel('EnvThreat');
    set(gca,'YTick',1:p.EnvThreat_Levels)
    set(gca, 'YTickLabel', ylabels)    
    xlabel('NTokens');
    set(gca,'XTick',1:p.NTokenPairs_Levels)
    set(gca, 'XTickLabel',xlabels)

end

