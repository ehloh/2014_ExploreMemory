function [trialpar norm p] = f_generate_taskstruc(p, col)
% [ trialpar norm p] = f_generate_taskstruc(p, col)
%
% Specify:
%     p          p.nreps_percell
%                 p.EnvThreat_Levels
%
%                 p.EnvThreat_pBomb (vector, optional)
%                     - vector specifying EnvThreat p(Bomb) probabilities
%                        (in ascending order of probability)
%                     - if empty, assume pBomb = increments of 1/p.EnvThreat_Levels,
%                         with maximal EnvThreat p(Bomb)=1. - i.e. (1:p.EnvThreat_Levels)/p.EnvThreat_Levels;
%                 p.NTokenPairs_Levels (optional)
%                     - if empty, assume same as EnvThreat
%                 p.NTokenPairs_Npairs (vector, optional)
%                     - if empty, assume 1:1: n_levels
%                 p.NSpaces_TokenPairs (optional)
%                     - if empty, assume same as max p.NTokenPairs_Npairs
%
%     col        col.EnvThreat
%                 col.NTokenPairs
%                 col.BombPresent
%                 col.BombActivated
%
%
%
% Output:
%
%     trialpar                    trial parameters (above cols filled)
%     norm.pActBomb       p(Activated Bomb) across design
%                                      - imagesc(norm.pActBomb) to view
%     norm.pBomb            p(Bomb) across design
%
%
% --------------------------------------------------------------------

% Execute to debug:  p.EnvThreat_Levels=2; p.nreps_percell=20; col.EnvThreat=2; col.NTokenPairs=3; col.BombPresent=4; col.BombActivated=5;  p.NTokenPairs_Npairs=[1 3];
% Execute to debug: p=p.design;


for o1=1:1 % Checks + Fill in the blanks for unspecified parameters
    if isfield(p,'EnvThreat_pBomb')==0
        p.EnvThreat_pBomb=(1:p.EnvThreat_Levels)/p.EnvThreat_Levels;
    end
    if isfield(p, 'NTokenPairs_Levels')==0
        p.NTokenPairs_Levels=p.EnvThreat_Levels;
    end
    if isfield(p, 'NTokenPairs_Npairs')==0
        p.NTokenPairs_Npairs=1:1:p.NTokenPairs_Levels;
    end
    if isfield(p,'NSpaces_TokenPairs')==0
        p.NSpaces_TokenPairs=p.NTokenPairs_Levels;
    end
    
    % Checks
    if length(p.EnvThreat_pBomb)~=p.EnvThreat_Levels;
        error('Invalid levels for NTokenPairs. length(p.EnvThreat_pBomb) must match p.EnvThreat_Levels');
    end
    if length(p.NTokenPairs_Npairs)~=p.NTokenPairs_Levels;
        error('Invalid levels for NTokenPairs. length(p.NTokenPairs_Npairs) must match p.NTokenPairs_Levels');
    end
    if length(p.NSpaces_TokenPairs)>p.NTokenPairs_Levels;
        error('Invalid NTokenPairs. p.NSpaces_TokenPairs cannot exceed p.NTokenPairs_Levels.');
    end
    if p.NSpaces_TokenPairs<p.NTokenPairs_Levels || p.NSpaces_TokenPairs<p.NTokenPairs_Npairs(end);
        error('Invalid p.NSpaces_TokenPairs - cannot be smaller than p.NTokenPairs_Levels  or the max p.NTokenPairs_Npairs');
    end
end

%% Assemble trial parameters

% Assemble Cell instructions: [1] EnvThreatLevel, [2] NTokenPairsLevel, [3] NTokenPairs [4] pBomb [5] pActBomb, [6] N Bombs, [7] N Activated Bombs
cells=zeros(p.EnvThreat_Levels*p.NTokenPairs_Levels,7);
cells(:,1)=ceil(1/p.EnvThreat_Levels:1/p.EnvThreat_Levels:p.EnvThreat_Levels)';
cells(:,2)=repmat((1:p.NTokenPairs_Levels)', p.EnvThreat_Levels,1);
for i=1:size(cells,1)
    cells(i,3)=p.NTokenPairs_Npairs(cells(i,2));
    cells(i,4)=p.EnvThreat_pBomb(cells(i,1));
    cells(i,5)=cells(i,4)*cells(i,3)/p.NSpaces_TokenPairs;
    cells(i,6)=round(p.nreps_percell*cells(i,4)); % round off
    cells(i,7)=round(p.nreps_percell*cells(i,5));
end

% Assemble trial parameters
trialpar=[];
for c=1:size(cells,1)
    wc.t=nan*zeros(p.nreps_percell, 4); % wc.t: [1] EnvThreat [2] NTokenPairs [3] Bomb present [4] Bomb activated
    wc.t(:,1)=cells(c,1);
    wc.t(:,2)=cells(c,3);
    wc.t(:,3)=0;
    wc.t(1:cells(c,6),3)=1;
    wc.t(:,4)=0;
    wc.t(1:cells(c,7),4)=1;
    
    % Write in requested format
    wc.tt(:,col.EnvThreat)=wc.t(:,1);
    wc.tt(:,col.NTokenPairs)=wc.t(:,2);
    wc.tt(:,col.BombPresent)=wc.t(:,3);
    wc.tt(:,col.BombActivated)=wc.t(:,4);
    
    trialpar=[trialpar; wc.tt];
    wc=[];
end

%% Assemble checks

% Assemble data-grid (for imagesc)
norm.pBomb=nan*zeros(p.EnvThreat_Levels, p.NTokenPairs_Levels);
norm.pActBomb=nan*zeros(p.EnvThreat_Levels, p.NTokenPairs_Levels);
for e=1:p.EnvThreat_Levels
    for n=1:p.NTokenPairs_Levels
        wg=trialpar(trialpar(:,col.EnvThreat)==e & trialpar(:,col.NTokenPairs)==p.NTokenPairs_Npairs(n), :);
        %
        norm.pBomb(p.EnvThreat_Levels+1-e,n)= sum(wg(:,col.BombPresent))/size(wg,1);
        norm.pActBomb(p.EnvThreat_Levels+1-e,n)= sum(wg(:,col.BombActivated))/size(wg,1);
    end
end

end

