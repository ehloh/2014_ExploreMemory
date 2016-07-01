function [model_defaults  par_specs models] = f_modelsettings(varargin) 
% [model_defaults  par_details models] = f_modelsettings(whichmodels, parsteps) 
% Generate model settings + fetch requested model details 'models'
% 
% Output variable 'models' lists details of each model -------------------------------------
%           Col 1: Model name
%           Col 2: No. free parameters
%           Col 3: Parameter names
%           Col 4: Constraint transform (constraint implemented in softmax/value fxns)
%           Col 5: Inverse constraint transform transform
%           Col 6: Parameter ranges
%
%   To run an ad-hoc model, add info to model_defaults (name, n free par, par
%   names), if there are new parameters, add into to par_transformations (par,
%   par_transformationsation in, par_transformationsation out)
%   (See function content/comment for details of model_defaults & par_transformations)
%
% --------------------------------------------------------------------------------------------

whichmodels=varargin{1};
if length(varargin)==2; parsteps=varargin{2};
    if parsteps==1; parsteps=2; end
else   parsteps=10; 
end

%% Parameter  specifications
%         Col 1: Parameter name
%         Col 2: Constraint transformations 
%         Col 3: Inverse constraint transformation
%         Col 4: Reasonable parameter range (in true parameter terms)
% 
% Constraint transforms (Col 2) are used to keep parameter values mathematically 
%   kosher, with respect to their usage in the softmax and value functions 
%   (e.g. probabilities between 0 and 1 etc). These transforms are applied within 
%   the softmax/value  functions, and listed here for reference. 
% 
% Inverse constraint transforms (Col 3) are applied if you want to use the existing
%   softmax/value function scripts to generate values/behaviour. Since the scripts
%   themselves apply the transformations, the 'fitted' parameter values must be 
%   inverse transformed first before they enter the softmax/value function scripts
%   themselves (or else, e.g. beta will be exponentiated twice).
% 
% Reasonable parameter range (Col 4), used to (a) seed model fitting iterations, 
%   and (b) specify parameter range in grid search. Because values are in'true'
%   parameter terms, values must be inverse constraint transformed before they
%   are fed (as seed values) into the softmax/value function scripts.
%   Ranges are not binding for the iterative model fits.
  
par_specs={       
     
    'b'     '20./(1+exp(-x))'        '- log( (20./x) - 1)'          [0.01    8];                    % b: softmax beta (beta>0)

    'p'     '1./(3+3.*exp(-x))'      '-log(  1./(3.*x)  -1 )'       [0.001    1/3-0.05];      % p: softmax epsilon ( 0 < epsilon < 1/3; sigmoid function)

    'm'     '20./(1+exp(-x))'        '- log( (20./x) - 1)'          [1/5        5] ;    % u: posterior probability power law distortion (0<m<20)
    
    'a'     '20./(1+exp(-x))'        '- log( (20./x) - 1)'   [0  20];   % a: post-Explore beta (0<a<20)
    
    'i'    '40./(1+exp(-x))  - 20'        '- log( (40./(x+20)) - 1)'   [-15 15] ;    % a: distort vAcceptGivenNoSee   (-20 < a < 20)
        
    'f'     'x'     'x'     [-20     -8];                                        % f: subjective value of loss
     
    'e'    'x'      'x'     [-2    4] ;                                       % e: exploration bonus, fixed magnitude

    'w'     'x' 'x'      [-1    1];                                 % w: exploration bonus, variable magnitude
                                                                                     % range =[0 1]. Range is  [0 5] for uncertainty multiplier (implemented below)
    
    'j'     '20./(1+exp(-x))'        '- log( (20./x) - 1)'          [1/5    5] ;    % j: envthreat probability power law distortion (0<j<20)
    
    
% % Redo ranges artificially  ###################################################

% % ####################################################################
    
    };

% The following are NOT free parameters; identity listed here for log.
%     See functions fpar_conflict and fpar_control for their specification
par_fixed={
        'v'     % VExplore:    'vw' = exploration bonus varies as a function of VExplore
        'u'     % Uncertainty/Entropy:  'uw' = exploration bonus varies as a function of Uncertainty
        'o'     % EntropyNTok:  'ow' = exploration bonus varies as a function of UncertaintyNTok
        'y'     % Variance in (binary) pActBomb:  'yw' = exploration bonus varies as a function of variance, i.e. pActBomb*(1-pActBomb)
        'k'     % Std dev (mean-variance theory) of outcome/EV:  'kw' = exploration bonus varies as a function of std dev. For details, see f_meanvar_quantities, fpar_conflict, fpar_control
        
        };
    
% Create parameter ranges from instructions
for p=1:size(par_specs,1)
   
    % Parameters that need to be evenly-spaced according to special scales (e.g. log scale)
    if strcmp(par_specs{p,1}, 'm')==1  % m = power law, even in log scales
        par_specs{p,4}= logspace( log10(par_specs{p,4}(1)) ,  log10(par_specs{p,4}(end)), parsteps);  % logspace is base 10. If steps are evenly spaced in 10, they will be evenly spaced in other bases as well
    else
        % Default parameter steps
          par_specs{p,4}= par_specs{p,4}(1)  : (par_specs{p,4}(end)-par_specs{p,4}(1))/(parsteps-1)  :    par_specs{p,4}(end);
    end
end


%% Model defaults: Free parameters 
%         Col 1: Model name
%         Col 2: # Free parameters
%         Col 3: Free parametes (in order, according to name)

model_defaults={
    'b01'          1       {'b'};                % Basic value function 
    'b02_f'       2        {'b'; 'f'};
    'b03_e'       2        {'b'; 'e'};
    'b04_uw'    2        {'b'; 'w'};
    'b05_vw'    2        {'b'; 'w'};
    'b06_ow'    2        {'b'; 'w'};
    'b07_fe'      3        {'b'; 'f'; 'e'};
    'b08_fuw'   3        {'b'; 'f'; 'w'};
    'b09_fvw'   3        {'b'; 'f'; 'w'};
    'b10_fow'   3        {'b'; 'f'; 'w'};
    'b11_euw'   3        {'b'; 'e'; 'w'};
    'b12_evw'   3        {'b'; 'e'; 'w'};
    'b13_eow'   3        {'b'; 'e'; 'w'};
    'b14_feuw'   4        {'b'; 'f'; 'e'; 'w'};
    'b15_fevw'   4        {'b'; 'f'; 'e'; 'w'};
    'b16_feow'   4        {'b'; 'f'; 'e'; 'w'};
    'b17_yw'    2        {'b'; 'w'};
    'b18_fyw'   3        {'b'; 'f'; 'w'};
    'b19_eyw'   3        {'b'; 'e'; 'w'};
    'b20_feyw'   4        {'b'; 'f'; 'e'; 'w'};
    'b21_kw'    2        {'b'; 'w'};
    'b22_fkw'   3        {'b'; 'f'; 'w'};
    'b23_ekw'   3        {'b'; 'e'; 'w'};
    'b24_fekw'   4        {'b'; 'f'; 'e'; 'w'};
    
    'bp01'          2       {'b'; 'p'};                % Epsilon (softmax lapse)
    'bp02_f'       3        {'b'; 'p'; 'f'};
    'bp03_e'       3        {'b'; 'p'; 'e'};
    'bp04_uw'    3        {'b'; 'p'; 'w'};
    'bp05_vw'    3        {'b'; 'p'; 'w'};
    'bp06_ow'    3        {'b'; 'p'; 'w'};
    'bp07_fe'      4        {'b'; 'p'; 'f'; 'e'};
    'bp08_fuw'   4        {'b'; 'p'; 'f'; 'w'};
    'bp09_fvw'   4        {'b'; 'p'; 'f'; 'w'};
    'bp10_fow'   4        {'b'; 'p'; 'f'; 'w'};
    'bp11_euw'   4        {'b'; 'p'; 'e'; 'w'};
    'bp12_evw'   4        {'b'; 'p'; 'e'; 'w'};
    'bp13_eow'   4        {'b'; 'p'; 'e'; 'w'};
    'bp14_feuw'   5        {'b'; 'p'; 'f'; 'e'; 'w'};
    'bp15_fevw'   5        {'b'; 'p'; 'f'; 'e'; 'w'};
    'bp16_feow'   5        {'b'; 'p'; 'f'; 'e'; 'w'};
    'bp17_yw'    3        {'b'; 'p'; 'w'};
    'bp18_fyw'   4        {'b'; 'p'; 'f'; 'w'};
    'bp19_eyw'   4        {'b'; 'p'; 'e'; 'w'};
    'bp20_feyw'   5        {'b'; 'p'; 'f'; 'e'; 'w'};
    'bp21_kw'    3        {'b'; 'p'; 'w'};
    'bp22_fkw'   4        {'b'; 'p'; 'f'; 'w'};
    'bp23_ekw'   4        {'b'; 'p'; 'e'; 'w'};
    'bp24_fekw'   5        {'b'; 'p'; 'f'; 'e'; 'w'};

    'bm01'         2       {'b'; 'm'};        % Unoptimal posterior probability distortion
    'bm02_f'      3        {'b'; 'm'; 'f'};
    'bm03_e'      3        {'b'; 'm'; 'e'};
    'bm04_uw'    3        {'b'; 'm'; 'w'};
    'bm05_vw'    3        {'b'; 'm'; 'w'};
    'bm06_ow'    3        {'b'; 'm'; 'w'};
    'bm07_fe'       4        {'b'; 'm'; 'f'; 'e'};
    'bm08_fuw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm09_fvw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm10_fow'    4        {'b'; 'm'; 'f'; 'w'};
    'bm11_euw'    4        {'b'; 'm'; 'e'; 'w'};
    'bm12_evw'    4        {'b'; 'm'; 'e'; 'w'};
    'bm13_eow'    4        {'b'; 'm'; 'e'; 'w'};
    'bm14_feuw'   5        {'b'; 'm'; 'f'; 'e'; 'w'};
    'bm15_fevw'   5        {'b'; 'm'; 'f'; 'e'; 'w'};
    'bm16_feow'   5        {'b'; 'm'; 'f'; 'e'; 'w'};
    'bm17_yw'    3        {'b'; 'm'; 'w'};
    'bm18_fyw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm19_eyw'    4        {'b'; 'm'; 'e'; 'w'};
    'bm20_feyw'   5        {'b'; 'm'; 'f'; 'e'; 'w'};    
    'bm21_kw'    3        {'b'; 'm'; 'w'};
    'bm22_fkw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm23_ekw'    4        {'b'; 'm'; 'e'; 'w'};
    'bm24_fekw'   5        {'b'; 'm'; 'f'; 'e'; 'w'};    

%     'ba01'         2       {'b'; 'a'};        % Distort vAcceeptGivenNoSee
%     'ba02_f'      3        {'b'; 'a'; 'f'};
%     'ba03_e'      3        {'b'; 'a'; 'e'};
%     'ba04_uw'    3        {'b'; 'a'; 'w'};
%     'ba05_vw'    3        {'b'; 'a'; 'w'};
%     'ba06_ow'    3        {'b'; 'a'; 'w'};
%     'ba07_fe'       4        {'b'; 'a'; 'f'; 'e'};
%     'ba08_fuw'    4        {'b'; 'a'; 'f'; 'w'};
%     'ba09_fvw'    4        {'b'; 'a'; 'f'; 'w'};
%     'ba10_fow'    4        {'b'; 'a'; 'f'; 'w'};
%     'ba11_euw'    4        {'b'; 'a'; 'e'; 'w'};
%     'ba12_evw'    4        {'b'; 'a'; 'e'; 'w'};
%     'ba13_eow'    4        {'b'; 'a'; 'e'; 'w'};
%     'ba14_feuw'   5        {'b'; 'a'; 'f'; 'e'; 'w'};
%     'ba15_fevw'   5        {'b'; 'a'; 'f'; 'e'; 'w'};
%     'ba16_feow'   5        {'b'; 'a'; 'f'; 'e'; 'w'};
%     'ba17_yw'    3        {'b'; 'a'; 'w'};
%     'ba18_fyw'    4        {'b'; 'a'; 'f'; 'w'};
%     'ba19_eyw'    4        {'b'; 'a'; 'e'; 'w'};
%     'ba20_feyw'   5        {'b'; 'a'; 'f'; 'e'; 'w'};
%     'ba21_kw'    3        {'b'; 'a'; 'w'};
%     'ba22_fkw'    4        {'b'; 'a'; 'f'; 'w'};
%     'ba23_ekw'    4        {'b'; 'a'; 'e'; 'w'};
%     'ba24_fekw'   5        {'b'; 'a'; 'f'; 'e'; 'w'};
    
    'bj01'         2       {'b'; 'j'};        % EnvThreat probability distortion
    'bj02_f'      3        {'b'; 'j'; 'f'};
    'bj03_e'      3        {'b'; 'j'; 'e'};   
    'bj04_uw'    3        {'b'; 'j'; 'w'};
    'bj05_vw'    3        {'b'; 'j'; 'w'};
    'bj06_ow'    3        {'b'; 'j'; 'w'};
    'bj07_fe'       4        {'b'; 'j'; 'f'; 'e'};
    'bj08_fuw'    4        {'b'; 'j'; 'f'; 'w'};
    'bj09_fvw'    4        {'b'; 'j'; 'f'; 'w'};
    'bj10_fow'    4        {'b'; 'j'; 'f'; 'w'};
    'bj11_euw'    4        {'b'; 'j'; 'e'; 'w'};
    'bj12_evw'    4        {'b'; 'j'; 'e'; 'w'};
    'bj13_eow'    4        {'b'; 'j'; 'e'; 'w'};
    'bj14_feuw'   5        {'b'; 'j'; 'f'; 'e'; 'w'};
    'bj15_fevw'   5        {'b'; 'j'; 'f'; 'e'; 'w'};
    'bj16_feow'   5        {'b'; 'j'; 'f'; 'e'; 'w'};
    'bj17_yw'    3        {'b'; 'j'; 'w'};
    'bj18_fyw'    4        {'b'; 'j'; 'f'; 'w'};
    'bj19_eyw'    4        {'b'; 'j'; 'e'; 'w'};
    'bj20_feyw'   5        {'b'; 'j'; 'f'; 'e'; 'w'};
    'bj21_kw'    3        {'b'; 'j'; 'w'};
    'bj22_fkw'    4        {'b'; 'j'; 'f'; 'w'};
    'bj23_ekw'    4        {'b'; 'j'; 'e'; 'w'};
    'bj24_fekw'   5        {'b'; 'j'; 'f'; 'e'; 'w'};
    
    % [3 base-param fams] #########################################
    
    'bpm01'         3       {'b'; 'p'; 'm'};        % Epsilon + Unoptimal posterior probability calculations (x Multiplier)
    'bpm02_f'      4        {'b'; 'p'; 'm'; 'f'};
    'bpm03_e'      4        {'b'; 'p'; 'm'; 'e'};
    'bpm04_uw'    4        {'b'; 'p'; 'm'; 'w'};
    'bpm05_vw'    4        {'b'; 'p'; 'm'; 'w'};
    'bpm06_ow'    4        {'b'; 'p'; 'm'; 'w'};
    'bpm07_fe'       5        {'b'; 'p'; 'm'; 'f'; 'e'};
    'bpm08_fuw'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm09_fvw'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm10_fow'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm11_euw'    5        {'b'; 'p'; 'm'; 'e'; 'w'};
    'bpm12_evw'    5        {'b'; 'p'; 'm'; 'e'; 'w'};
    'bpm13_eow'    5        {'b'; 'p'; 'm'; 'e'; 'w'};
    'bpm14_feuw'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    'bpm15_fevw'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    'bpm16_feow'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    'bpm17_yw'    4        {'b'; 'p'; 'm'; 'w'};
    'bpm18_fyw'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm19_eyw'    5        {'b'; 'p'; 'm'; 'e'; 'w'};
    'bpm20_feyw'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    'bpm21_kw'    4        {'b'; 'p'; 'm'; 'w'};
    'bpm22_fkw'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm23_ekw'    5        {'b'; 'p'; 'm'; 'e'; 'w'};
    'bpm24_fekw'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'}; 
    %
%     'bma01'         3       {'b'; 'm'; 'a'};        % Unoptimal posterior probability + vAccGivNoSee bonus
%     'bma02_f'      4        {'b'; 'm'; 'a'; 'f'};
%     'bma03_e'      4        {'b'; 'm'; 'a'; 'e'};
%     'bma04_uw'    4        {'b'; 'm'; 'a'; 'w'};
%     'bma05_vw'    4        {'b'; 'm'; 'a'; 'w'};
%     'bma06_ow'    4        {'b'; 'm'; 'a'; 'w'};
%     'bma07_fe'       5        {'b'; 'm'; 'a'; 'f'; 'e'};
%     'bma08_fuw'    5        {'b'; 'm'; 'a'; 'f'; 'w'};
%     'bma09_fvw'    5        {'b'; 'm'; 'a'; 'f'; 'w'};
%     'bma10_fow'    5        {'b'; 'm'; 'a'; 'f'; 'w'};
%     'bma11_euw'    5        {'b'; 'm'; 'a'; 'e'; 'w'};
%     'bma12_evw'    5        {'b'; 'm'; 'a'; 'e'; 'w'};
%     'bma13_eow'    5        {'b'; 'm'; 'a'; 'e'; 'w'};
%     'bma14_feuw'   6        {'b'; 'm'; 'a'; 'f'; 'e'; 'w'};
%     'bma15_fevw'   6        {'b'; 'm'; 'a'; 'f'; 'e'; 'w'};
%     'bma16_feow'   6        {'b'; 'm'; 'a'; 'f'; 'e'; 'w'};
%     'bma17_yw'    4        {'b'; 'm'; 'a'; 'w'};
%     'bma18_fyw'    5        {'b'; 'm'; 'a'; 'f'; 'w'};
%     'bma19_eyw'    5        {'b'; 'm'; 'a'; 'e'; 'w'};
%     'bma20_feyw'   6        {'b'; 'm'; 'a'; 'f'; 'e'; 'w'};
%     'bma21_kw'    4        {'b'; 'm'; 'a'; 'w'};
%     'bma22_fkw'    5        {'b'; 'm'; 'a'; 'f'; 'w'};
%     'bma23_ekw'    5        {'b'; 'm'; 'a'; 'e'; 'w'};
%     'bma24_fekw'   6        {'b'; 'm'; 'a'; 'f'; 'e'; 'w'};
    
    %
    'bjm01'         3       {'b'; 'j'; 'm'};        % Epsilon + EnvThreat distortion + posterior prob distortion
    'bjm02_f'      4        {'b'; 'j'; 'm'; 'f'};
    'bjm03_e'      4        {'b'; 'j'; 'm'; 'e'};
    'bjm04_uw'    4        {'b'; 'j'; 'm'; 'w'};
    'bjm05_vw'    4        {'b'; 'j'; 'm'; 'w'};
    'bjm06_ow'    4        {'b'; 'j'; 'm'; 'w'};
    'bjm07_fe'       5        {'b'; 'j'; 'm'; 'f'; 'e'};
    'bjm08_fuw'    5        {'b'; 'j'; 'm'; 'f'; 'w'};
    'bjm09_fvw'    5        {'b'; 'j'; 'm'; 'f'; 'w'};
    'bjm10_fow'    5        {'b'; 'j'; 'm'; 'f'; 'w'};
    'bjm11_euw'    5        {'b'; 'j'; 'm'; 'e'; 'w'};
    'bjm12_evw'    5        {'b'; 'j'; 'm'; 'e'; 'w'};
    'bjm13_eow'    5        {'b'; 'j'; 'm'; 'e'; 'w'};
    'bjm14_feuw'   6        {'b'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bjm15_fevw'   6        {'b'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bjm16_feow'   6        {'b'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bjm17_yw'    4        {'b'; 'j'; 'm'; 'w'};
    'bjm18_fyw'    5        {'b'; 'j'; 'm'; 'f'; 'w'};
    'bjm19_eyw'    5        {'b'; 'j'; 'm'; 'e'; 'w'};
    'bjm20_feyw'   6        {'b'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bjm21_kw'    4        {'b'; 'j'; 'm'; 'w'};
    'bjm22_fkw'    5        {'b'; 'j'; 'm'; 'f'; 'w'};
    'bjm23_ekw'    5        {'b'; 'j'; 'm'; 'e'; 'w'};
    'bjm24_fekw'   6        {'b'; 'j'; 'm'; 'f'; 'e'; 'w'};
    %
    'bpj01'         3       {'b'; 'p'; 'j'};        % Epsilon + Unoptijal posterior probability calculations (x jultiplier)
    'bpj02_f'      4        {'b'; 'p'; 'j'; 'f'};
    'bpj03_e'      4        {'b'; 'p'; 'j'; 'e'};
    'bpj04_uw'    4        {'b'; 'p'; 'j'; 'w'};
    'bpj05_vw'    4        {'b'; 'p'; 'j'; 'w'};
    'bpj06_ow'    4        {'b'; 'p'; 'j'; 'w'};
    'bpj07_fe'       5        {'b'; 'p'; 'j'; 'f'; 'e'};
    'bpj08_fuw'    5        {'b'; 'p'; 'j'; 'f'; 'w'};
    'bpj09_fvw'    5        {'b'; 'p'; 'j'; 'f'; 'w'};
    'bpj10_fow'    5        {'b'; 'p'; 'j'; 'f'; 'w'};
    'bpj11_euw'    5        {'b'; 'p'; 'j'; 'e'; 'w'};
    'bpj12_evw'    5        {'b'; 'p'; 'j'; 'e'; 'w'};
    'bpj13_eow'    5        {'b'; 'p'; 'j'; 'e'; 'w'};
    'bpj14_feuw'   6        {'b'; 'p'; 'j'; 'f'; 'e'; 'w'};
    'bpj15_fevw'   6        {'b'; 'p'; 'j'; 'f'; 'e'; 'w'};
    'bpj16_feow'   6        {'b'; 'p'; 'j'; 'f'; 'e'; 'w'};
    'bpj17_yw'    4        {'b'; 'p'; 'j'; 'w'};
    'bpj18_fyw'    5        {'b'; 'p'; 'j'; 'f'; 'w'};
    'bpj19_eyw'    5        {'b'; 'p'; 'j'; 'e'; 'w'};
    'bpj20_feyw'   6        {'b'; 'p'; 'j'; 'f'; 'e'; 'w'};
    'bpj21_kw'    4        {'b'; 'p'; 'j'; 'w'};
    'bpj22_fkw'    5        {'b'; 'p'; 'j'; 'f'; 'w'};
    'bpj23_ekw'    5        {'b'; 'p'; 'j'; 'e'; 'w'};
    'bpj24_fekw'   6        {'b'; 'p'; 'j'; 'f'; 'e'; 'w'};
    %
%     'bji01'         3       {'b'; 'j'; 'i'};        % Epsilon + EnvThreat distortion + posterior prob distortion
%     'bji02_f'      4        {'b'; 'j'; 'i'; 'f'};
%     'bji03_e'      4        {'b'; 'j'; 'i'; 'e'};
%     'bji04_uw'    4        {'b'; 'j'; 'i'; 'w'};
%     'bji05_vw'    4        {'b'; 'j'; 'i'; 'w'};
%     'bji06_ow'    4        {'b'; 'j'; 'i'; 'w'};
%     'bji07_fe'       5        {'b'; 'j'; 'i'; 'f'; 'e'};
%     'bji08_fuw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
%     'bji09_fvw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
%     'bji10_fow'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
%     'bji11_euw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
%     'bji12_evw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
%     'bji13_eow'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
%     'bji14_feuw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
%     'bji15_fevw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
%     'bji16_feow'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
%     'bji17_yw'    4        {'b'; 'j'; 'i'; 'w'};
%     'bji18_fyw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
%     'bji19_eyw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
%     'bji20_feyw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
%     'bji21_kw'    4        {'b'; 'j'; 'i'; 'w'};
%     'bji22_fkw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
%     'bji23_ekw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
%     'bji24_fekw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
%     
    %  #########################################
    %
    'bpjm01'         4       {'b'; 'p'; 'j';'m'};        % Epsilon + Unoptijal posterior probability calculations (x jultiplier)
    'bpjm02_f'      5        {'b'; 'p'; 'j'; 'm'; 'f'};
    'bpjm03_e'      5        {'b'; 'p'; 'j'; 'm'; 'e'};
    'bpjm04_uw'    5        {'b'; 'p'; 'j'; 'm'; 'w'};
    'bpjm05_vw'    5        {'b'; 'p'; 'j'; 'm'; 'w'};
    'bpjm06_ow'    5        {'b'; 'p'; 'j'; 'm'; 'w'};
    'bpjm07_fe'       6        {'b'; 'p'; 'j'; 'm'; 'f'; 'e'};
    'bpjm08_fuw'    6        {'b'; 'p'; 'j'; 'm'; 'f'; 'w'};
    'bpjm09_fvw'    6        {'b'; 'p'; 'j'; 'm'; 'f'; 'w'};
    'bpjm10_fow'    6        {'b'; 'p'; 'j'; 'm'; 'f'; 'w'};
    'bpjm11_euw'    6        {'b'; 'p'; 'j'; 'm'; 'e'; 'w'};
    'bpjm12_evw'    6        {'b'; 'p'; 'j'; 'm'; 'e'; 'w'};
    'bpjm13_eow'    6        {'b'; 'p'; 'j'; 'm'; 'e'; 'w'};
    'bpjm14_feuw'   7        {'b'; 'p'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bpjm15_fevw'   7        {'b'; 'p'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bpjm16_feow'   7        {'b'; 'p'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bpjm17_yw'    5        {'b'; 'p'; 'j'; 'm'; 'w'};
    'bpjm18_fyw'    6        {'b'; 'p'; 'j'; 'm'; 'f'; 'w'};
    'bpjm19_eyw'    6        {'b'; 'p'; 'j'; 'm'; 'e'; 'w'};
    'bpjm20_feyw'   7        {'b'; 'p'; 'j'; 'm'; 'f'; 'e'; 'w'};
    'bpjm21_kw'    5        {'b'; 'p'; 'j'; 'm'; 'w'};
    'bpjm22_fkw'    6        {'b'; 'p'; 'j'; 'm'; 'f'; 'w'};
    'bpjm23_ekw'    6        {'b'; 'p'; 'j'; 'm'; 'e'; 'w'};
    'bpjm24_fekw'   7        {'b'; 'p'; 'j'; 'm'; 'f'; 'e'; 'w'};
    
    % ##########################################
    % ##########################################
    % ##########################################
    
    'bi01'         2       {'b'; 'i'};        % 
    'bi02_f'      3        {'b'; 'i'; 'f'};
    'bi03_e'      3        {'b'; 'i'; 'e'};
    'bi04_uw'    3        {'b'; 'i'; 'w'};
    'bi05_vw'    3        {'b'; 'i'; 'w'};
    'bi06_ow'    3        {'b'; 'i'; 'w'};
    'bi07_fe'       4        {'b'; 'i'; 'f'; 'e'};
    'bi08_fuw'    4        {'b'; 'i'; 'f'; 'w'};
    'bi09_fvw'    4        {'b'; 'i'; 'f'; 'w'};
    'bi10_fow'    4        {'b'; 'i'; 'f'; 'w'};
    'bi11_euw'    4        {'b'; 'i'; 'e'; 'w'};
    'bi12_evw'    4        {'b'; 'i'; 'e'; 'w'};
    'bi13_eow'    4        {'b'; 'i'; 'e'; 'w'};
    'bi14_feuw'   5        {'b'; 'i'; 'f'; 'e'; 'w'};
    'bi15_fevw'   5        {'b'; 'i'; 'f'; 'e'; 'w'};
    'bi16_feow'   5        {'b'; 'i'; 'f'; 'e'; 'w'};
    'bi17_yw'    3        {'b'; 'i'; 'w'};
    'bi18_fyw'    4        {'b'; 'i'; 'f'; 'w'};
    'bi19_eyw'    4        {'b'; 'i'; 'e'; 'w'};
    'bi20_feyw'   5        {'b'; 'i'; 'f'; 'e'; 'w'};
    'bi21_kw'    3        {'b'; 'i'; 'w'};
    'bi22_fkw'    4        {'b'; 'i'; 'f'; 'w'};
    'bi23_ekw'    4        {'b'; 'i'; 'e'; 'w'};
    'bi24_fekw'   5        {'b'; 'i'; 'f'; 'e'; 'w'};
    
    %     -------------------------------------------------------------------------------------
    
    'bpi01'         3       {'b'; 'p'; 'i'};        % Epsilon + viccGivNoSee bonus
    'bpi02_f'      4        {'b'; 'p'; 'i'; 'f'};
    'bpi03_e'      4        {'b'; 'p'; 'i'; 'e'};
    'bpi04_uw'    4        {'b'; 'p'; 'i'; 'w'};
    'bpi05_vw'    4        {'b'; 'p'; 'i'; 'w'};
    'bpi06_ow'    4        {'b'; 'p'; 'i'; 'w'};
    'bpi07_fe'       5        {'b'; 'p'; 'i'; 'f'; 'e'};
    'bpi08_fuw'    5        {'b'; 'p'; 'i'; 'f'; 'w'};
    'bpi09_fvw'    5        {'b'; 'p'; 'i'; 'f'; 'w'};
    'bpi10_fow'    5        {'b'; 'p'; 'i'; 'f'; 'w'};
    'bpi11_euw'    5        {'b'; 'p'; 'i'; 'e'; 'w'};
    'bpi12_evw'    5        {'b'; 'p'; 'i'; 'e'; 'w'};
    'bpi13_eow'    5        {'b'; 'p'; 'i'; 'e'; 'w'};
    'bpi14_feuw'   6        {'b'; 'p'; 'i'; 'f'; 'e'; 'w'};
    'bpi15_fevw'   6        {'b'; 'p'; 'i'; 'f'; 'e'; 'w'};
    'bpi16_feow'   6        {'b'; 'p'; 'i'; 'f'; 'e'; 'w'};    
    'bpi17_yw'    4        {'b'; 'p'; 'i'; 'w'};
    'bpi18_fyw'    5        {'b'; 'p'; 'i'; 'f'; 'w'};
    'bpi19_eyw'    5        {'b'; 'p'; 'i'; 'e'; 'w'};
    'bpi20_feyw'   6        {'b'; 'p'; 'i'; 'f'; 'e'; 'w'};
    'bpi21_kw'    4        {'b'; 'p'; 'i'; 'w'};
    'bpi22_fkw'    5        {'b'; 'p'; 'i'; 'f'; 'w'};
    'bpi23_ekw'    5        {'b'; 'p'; 'i'; 'e'; 'w'};
    'bpi24_fekw'   6        {'b'; 'p'; 'i'; 'f'; 'e'; 'w'};
    
    'bmi01'         3       {'b'; 'm'; 'i'};        % Unoptimil posterior probibility + viccGivNoSee bonus
    'bmi02_f'      4        {'b'; 'm'; 'i'; 'f'};
    'bmi03_e'      4        {'b'; 'm'; 'i'; 'e'};
    'bmi04_uw'    4        {'b'; 'm'; 'i'; 'w'};
    'bmi05_vw'    4        {'b'; 'm'; 'i'; 'w'};
    'bmi06_ow'    4        {'b'; 'm'; 'i'; 'w'};
    'bmi07_fe'       5        {'b'; 'm'; 'i'; 'f'; 'e'};
    'bmi08_fuw'    5        {'b'; 'm'; 'i'; 'f'; 'w'};
    'bmi09_fvw'    5        {'b'; 'm'; 'i'; 'f'; 'w'};
    'bmi10_fow'    5        {'b'; 'm'; 'i'; 'f'; 'w'};
    'bmi11_euw'    5        {'b'; 'm'; 'i'; 'e'; 'w'};
    'bmi12_evw'    5        {'b'; 'm'; 'i'; 'e'; 'w'};
    'bmi13_eow'    5        {'b'; 'm'; 'i'; 'e'; 'w'};
    'bmi14_feuw'   6        {'b'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bmi15_fevw'   6        {'b'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bmi16_feow'   6        {'b'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bmi17_yw'    4        {'b'; 'm'; 'i'; 'w'};
    'bmi18_fyw'    5        {'b'; 'm'; 'i'; 'f'; 'w'};
    'bmi19_eyw'    5        {'b'; 'm'; 'i'; 'e'; 'w'};
    'bmi20_feyw'   6        {'b'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bmi21_kw'    4        {'b'; 'm'; 'i'; 'w'};
    'bmi22_fkw'    5        {'b'; 'm'; 'i'; 'f'; 'w'};
    'bmi23_ekw'    5        {'b'; 'm'; 'i'; 'e'; 'w'};
    'bmi24_fekw'   6        {'b'; 'm'; 'i'; 'f'; 'e'; 'w'};
    %
    'bpmi01'         4       {'b'; 'p'; 'm'; 'i'};        % Epsilon + Unoptimil posterior probibility cilculition + viccGivNoSee bonus
    'bpmi02_f'      5        {'b'; 'p'; 'm'; 'i'; 'f'};
    'bpmi03_e'      5        {'b'; 'p'; 'm'; 'i'; 'e'};
    'bpmi04_uw'    5        {'b'; 'p'; 'm'; 'i'; 'w'};
    'bpmi05_vw'    5        {'b'; 'p'; 'm'; 'i'; 'w'};
    'bpmi06_ow'    5        {'b'; 'p'; 'm'; 'i'; 'w'};
    'bpmi07_fe'       6        {'b'; 'p'; 'm'; 'i'; 'f'; 'e'};
    'bpmi08_fuw'    6        {'b'; 'p'; 'm'; 'i'; 'f'; 'w'};
    'bpmi09_fvw'    6        {'b'; 'p'; 'm'; 'i'; 'f'; 'w'};
    'bpmi10_fow'    6        {'b'; 'p'; 'm'; 'i'; 'f'; 'w'};
    'bpmi11_euw'    6        {'b'; 'p'; 'm'; 'i'; 'e'; 'w'};
    'bpmi12_evw'    6        {'b'; 'p'; 'm'; 'i'; 'e'; 'w'};
    'bpmi13_eow'    6        {'b'; 'p'; 'm'; 'i'; 'e'; 'w'};
    'bpmi14_feuw'   7        {'b'; 'p'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpmi15_fevw'   7        {'b'; 'p'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpmi16_feow'   7        {'b'; 'p'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpmi17_yw'    5        {'b'; 'p'; 'm'; 'i'; 'w'};
    'bpmi18_fyw'    6        {'b'; 'p'; 'm'; 'i'; 'f'; 'w'};
    'bpmi19_eyw'    6        {'b'; 'p'; 'm'; 'i'; 'e'; 'w'};
    'bpmi20_feyw'   7        {'b'; 'p'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpmi21_kw'    5        {'b'; 'p'; 'm'; 'i'; 'w'};
    'bpmi22_fkw'    6        {'b'; 'p'; 'm'; 'i'; 'f'; 'w'};
    'bpmi23_ekw'    6        {'b'; 'p'; 'm'; 'i'; 'e'; 'w'};
    'bpmi24_fekw'   7        {'b'; 'p'; 'm'; 'i'; 'f'; 'e'; 'w'};
    
    % -------
    
    'bji01'         3       {'b'; 'j'; 'i'};        % j i
    'bji02_f'      4        {'b'; 'j'; 'i'; 'f'};
    'bji03_e'      4        {'b'; 'j'; 'i'; 'e'};
    'bji04_uw'    4        {'b'; 'j'; 'i'; 'w'};
    'bji05_vw'    4        {'b'; 'j'; 'i'; 'w'};
    'bji06_ow'    4        {'b'; 'j'; 'i'; 'w'};
    'bji07_fe'       5        {'b'; 'j'; 'i'; 'f'; 'e'};
    'bji08_fuw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
    'bji09_fvw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
    'bji10_fow'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
    'bji11_euw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
    'bji12_evw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
    'bji13_eow'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
    'bji14_feuw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bji15_fevw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bji16_feow'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bji17_yw'    4        {'b'; 'j'; 'i'; 'w'};
    'bji18_fyw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
    'bji19_eyw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
    'bji20_feyw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bji21_kw'    4        {'b'; 'j'; 'i'; 'w'};
    'bji22_fkw'    5        {'b'; 'j'; 'i'; 'f'; 'w'};
    'bji23_ekw'    5        {'b'; 'j'; 'i'; 'e'; 'w'};
    'bji24_fekw'   6        {'b'; 'j'; 'i'; 'f'; 'e'; 'w'};
    
    'bpji01'         4       {'b'; 'p'; 'j'; 'i'};        % p j i
    'bpji02_f'      5        {'b'; 'p'; 'j'; 'i'; 'f'};
    'bpji03_e'      5        {'b'; 'p'; 'j'; 'i'; 'e'};
    'bpji04_uw'    5        {'b'; 'p'; 'j'; 'i'; 'w'};
    'bpji05_vw'    5        {'b'; 'p'; 'j'; 'i'; 'w'};
    'bpji06_ow'    5        {'b'; 'p'; 'j'; 'i'; 'w'};
    'bpji07_fe'       6        {'b'; 'p'; 'j'; 'i'; 'f'; 'e'};
    'bpji08_fuw'    6        {'b'; 'p'; 'j'; 'i'; 'f'; 'w'};
    'bpji09_fvw'    6        {'b'; 'p'; 'j'; 'i'; 'f'; 'w'};
    'bpji10_fow'    6        {'b'; 'p'; 'j'; 'i'; 'f'; 'w'};
    'bpji11_euw'    6        {'b'; 'p'; 'j'; 'i'; 'e'; 'w'};
    'bpji12_evw'    6        {'b'; 'p'; 'j'; 'i'; 'e'; 'w'};
    'bpji13_eow'    6        {'b'; 'p'; 'j'; 'i'; 'e'; 'w'};
    'bpji14_feuw'   7        {'b'; 'p'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bpji15_fevw'   7        {'b'; 'p'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bpji16_feow'   7        {'b'; 'p'; 'j'; 'i'; 'f'; 'e'; 'w'};    
    'bpji17_yw'    5        {'b'; 'p'; 'j'; 'i'; 'w'};
    'bpji18_fyw'    6        {'b'; 'p'; 'j'; 'i'; 'f'; 'w'};
    'bpji19_eyw'    6        {'b'; 'p'; 'j'; 'i'; 'e'; 'w'};
    'bpji20_feyw'   7        {'b'; 'p'; 'j'; 'i'; 'f'; 'e'; 'w'};
    'bpji21_kw'    5        {'b'; 'p'; 'j'; 'i'; 'w'};
    'bpji22_fkw'    6        {'b'; 'p'; 'j'; 'i'; 'f'; 'w'};
    'bpji23_ekw'    6        {'b'; 'p'; 'j'; 'i'; 'e'; 'w'};
    'bpji24_fekw'   7        {'b'; 'p'; 'j'; 'i'; 'f'; 'e'; 'w'};
    
    'bjmi01'         4       {'b'; 'j'; 'm'; 'i'};      
    'bjmi02_f'      5        {'b'; 'j'; 'm'; 'i'; 'f'};
    'bjmi03_e'      5        {'b'; 'j'; 'm'; 'i'; 'e'};
    'bjmi04_uw'    5        {'b'; 'j'; 'm'; 'i'; 'w'};
    'bjmi05_vw'    5        {'b'; 'j'; 'm'; 'i'; 'w'};
    'bjmi06_ow'    5        {'b'; 'j'; 'm'; 'i'; 'w'};
    'bjmi07_fe'       6        {'b'; 'j'; 'm'; 'i'; 'f'; 'e'};
    'bjmi08_fuw'    6        {'b'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bjmi09_fvw'    6        {'b'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bjmi10_fow'    6        {'b'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bjmi11_euw'    6        {'b'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bjmi12_evw'    6        {'b'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bjmi13_eow'    6        {'b'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bjmi14_feuw'   7        {'b'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bjmi15_fevw'   7        {'b'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bjmi16_feow'   7        {'b'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bjmi17_yw'    5        {'b'; 'j'; 'm'; 'i'; 'w'};
    'bjmi18_fyw'    6        {'b'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bjmi19_eyw'    6        {'b'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bjmi20_feyw'   7        {'b'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bjmi21_kw'    5        {'b'; 'j'; 'm'; 'i'; 'w'};
    'bjmi22_fkw'    6        {'b'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bjmi23_ekw'    6        {'b'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bjmi24_fekw'   7        {'b'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    %
    'bpjmi01'         5       {'b'; 'p'; 'j'; 'm'; 'i'};        
    'bpjmi02_f'      6        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'};
    'bpjmi03_e'      6        {'b'; 'p'; 'j'; 'm'; 'i'; 'e'};
    'bpjmi04_uw'    6        {'b'; 'p'; 'j'; 'm'; 'i'; 'w'};
    'bpjmi05_vw'    6        {'b'; 'p'; 'j'; 'm'; 'i'; 'w'};
    'bpjmi06_ow'    6        {'b'; 'p'; 'j'; 'm'; 'i'; 'w'};
    'bpjmi07_fe'       7        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'e'};
    'bpjmi08_fuw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bpjmi09_fvw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bpjmi10_fow'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bpjmi11_euw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bpjmi12_evw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bpjmi13_eow'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bpjmi14_feuw'   8        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpjmi15_fevw'   8        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpjmi16_feow'   8        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpjmi17_yw'    6        {'b'; 'p'; 'j'; 'm'; 'i'; 'w'};
    'bpjmi18_fyw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bpjmi19_eyw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bpjmi20_feyw'   8        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    'bpjmi21_kw'    6        {'b'; 'p'; 'j'; 'm'; 'i'; 'w'};
    'bpjmi22_fkw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'w'};
    'bpjmi23_ekw'    7        {'b'; 'p'; 'j'; 'm'; 'i'; 'e'; 'w'};
    'bpjmi24_fekw'   8        {'b'; 'p'; 'j'; 'm'; 'i'; 'f'; 'e'; 'w'};
    
    
    % -----------------------------------------------
    'allpars'            7        {'b'; 'p'; 'j'; 'm';'a'; 'f'; 'e'; 'w'; 'i'};
    
    };

% Correct no. of params?
for m=1:size(model_defaults,1)
    model_defaults{m,2}=length(model_defaults{m,3});
end

%% Generate details for requested models: 'models'
%       See function description (above) for description of 'models'
% keyboard

models=cell(length(whichmodels),6); % Read details
for m=1:length(whichmodels)
    if sum(strcmp(model_defaults(:,1), whichmodels{m}))~=0
        wm.modnum=find(strcmp(model_defaults(:,1), whichmodels{m}));
        models{m,1}=model_defaults{wm.modnum,1};
        models{m,2}=model_defaults{wm.modnum,2};
        models{m,3}=model_defaults{wm.modnum,3};
        
        % Fetch details for individual parameters
        for p=1:model_defaults{wm.modnum,2}
            wm.parnum=find(strcmp(par_specs(:,1), model_defaults{wm.modnum,3}{p}));
            models{m,4}{p}=par_specs{wm.parnum,2};
            models{m,5}{p}=par_specs{wm.parnum,3};
            models{m,6}{p}=par_specs{wm.parnum,4};
            
            % Uncertainty-multiplied exploration bonus: scale
            if isempty(strfind(models{m,1}, 'uw'))==0 & strcmp(models{m,3}{p}, 'w')
                models{m,6}{p}=models{m,6}{p}*5;
            end
        end
    else error(['Cannot find requested model in f_modelsettings: ' whichmodels{m}])
    end
end


end

