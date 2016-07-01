function [ pprob] = f_posteriorprob( et, ntok)
%  p(ActBomb| ~See) =  (E*N)/ (2B - N*E)
%   where E= EnvThreat, N= NTok, B= Base # tokens, p(See|ActBomb)=0.5

pprob=(et.*ntok) ./  (2*8-(et.*ntok));

end