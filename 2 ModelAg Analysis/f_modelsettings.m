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
if length(varargin)==2
    parsteps=varargin{2};
else
    parsteps=10; % disp('No. of parameter steps not specified. Default=10.')
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

    'bm01'         2       {'b'; 'm'};        % Unoptimal posterior probability calculations (x Multiplier)
    'bm02_f'      3        {'b'; 'm'; 'f'};
    'bm03_e'      3        {'b'; 'm'; 'e'};
    'bm04_uw'    3        {'b'; 'm'; 'w'};
    'bm05_vw'    3        {'b'; 'm'; 'w'};
    'bm06_ow'    3        {'b'; 'm'; 'w'};
    'bm07_fe'       4        {'b'; 'm'; 'f'; 'e'};
    'bm08_fuw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm09_fvw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm10_fow'    4        {'b'; 'm'; 'f'; 'w'};
    'bm11_euw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm12_evw'    4        {'b'; 'm'; 'f'; 'w'};
    'bm13_eow'    4        {'b'; 'm'; 'f'; 'w'};
    'bm14_feuw'   5        {'b'; 'm'; 'f'; 'e'; 'w'};
    'bm15_fevw'   5        {'b'; 'm'; 'f'; 'e'; 'w'};
    'bm16_feow'   5        {'b'; 'm'; 'f'; 'e'; 'w'};
    
    
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
%     
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
    'bpm11_euw'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm12_evw'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm13_eow'    5        {'b'; 'p'; 'm'; 'f'; 'w'};
    'bpm14_feuw'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    'bpm15_fevw'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    'bpm16_feow'   6        {'b'; 'p'; 'm'; 'f'; 'e'; 'w'};
    
    % ECONOMIC MODELS -----------------------------------------------
    
    'batl01'         4       {'b'; 'a'; 't'; 'l'};        %  Economic alpha (gains power distortion) and lambda
    'batl02_f'      5        {'b'; 'a'; 't'; 'l'; 'f'};
    'batl03_e'      5        {'b'; 'a'; 't'; 'l'; 'e'};
    'batl04_uw'    5        {'b'; 'a'; 't'; 'l'; 'w'};
    'batl05_vw'    5        {'b'; 'a'; 't'; 'l'; 'w'};
    'batl06_ow'    5        {'b'; 'a'; 't'; 'l'; 'w'};
    'batl07_fe'       6        {'b'; 'a'; 't'; 'l'; 'f'; 'e'};
    'batl08_fuw'    6        {'b'; 'a'; 't'; 'l'; 'f'; 'w'};
    'batl09_fvw'    6        {'b'; 'a'; 't'; 'l'; 'f'; 'w'};
    'batl10_fow'    6        {'b'; 'a'; 't'; 'l'; 'f'; 'w'};
    'batl11_euw'    6        {'b'; 'a'; 't'; 'l'; 'f'; 'w'};
    'batl12_evw'    6        {'b'; 'a'; 't'; 'l'; 'f'; 'w'};
    'batl13_eow'    6        {'b'; 'a'; 't'; 'l'; 'f'; 'w'};
    'batl14_feuw'   7        {'b'; 'a'; 't'; 'l'; 'f'; 'e'; 'w'};
    'batl15_fevw'   7        {'b'; 'a'; 't'; 'l'; 'f'; 'e'; 'w'};
    'batl16_feow'   7        {'b'; 'a'; 't'; 'l'; 'f'; 'e'; 'w'};
    
    % -----------------------------------------------
    'allpars'            7        {'b'; 'p'; 'm'; 'j'; 'f'; 'e'; 'w'};
    
    };

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
    
    'b'     'exp(x)'      'log(x)'      [0.01    6];                                             % b: softmax beta (beta>0)

     'p'     '1./(3+3.*exp(-x))'      '-log(  1./(3.*x)  -1 )'       [0.001    0.005];             % p: softmax epsilon ( 0 < epsilon < 1/3; sigmoid function)

    'm'     'exp(x)'     'log(x)'        [1/8        5] ;                                        % u: posterior probability power law distortion
    
    'f'     'x'     'x'     [-20     -8];                                        % f: subjective value of loss
    
    'e'    'x'      'x'     [-13    2] ;                                       % e: exploration bonus, fixed magnitude
                                                                                    %
    
    'w'     'x' 'x'      [0.5    2.8];                                 % w: exploration bonus, variable magnitude
                                                                          %             range =[0 1]. Range is  [0 5] for uncertainty multiplier (implemented below)
	'j'     '20./(1+exp(-x))'        '- log( (20./x) - 1)'          [1/5    5] ;    % j: envthreat probability power law distortion (0<j<20)

    % ----------------------------
    
    'a'    '1/(1+exp(-x))'      '-log((1-x)/x)'      [0.001    1];   % a: alpha, power for gain; 0<a<1
    
    't'    '1/(1+exp(-x))'      '-log((1-x)/x)'      [0.001    1];   % t: beta, power for losses; 0<t<1
    
    'l'     'exp(x)+1'      'log(x-1)'        [1    10];    % l: lambda, loss-aversion parameter; lambda>1
    
    
    
    };


% The following are NOT free parameters; identity listed here for log.
%     See functions fpar_conflict and fpar_control for their specification
par_fixed={
    
        'v'     % VExplore:    'vw' = exploration bonus varies as a function of VExplore
        'u'     % Uncertainty/Entropy:  'uw' = exploration bonus varies as a function of Uncertainty
        'o'     % EntropyNTok:  'ow' = exploration bonus varies as a function of Uncertainty

        };
    
% Create parameter ranges from instructions
for p=1:size(par_specs,1)
   
    % Parameters that need to be evenly-spaced according to special scales (e.g. log scale)
    if strcmp(par_specs{p,1}, 'm')==1;  % m = power law, even in log scales
        par_specs{p,4}= logspace( log10(par_specs{p,4}(1)) ,  log10(par_specs{p,4}(end)), parsteps);  % logspace is base 10. If steps are evenly spaced in 10, they will be evenly spaced in other bases as well
    else
        
        % Default parameter steps
          par_specs{p,4}= par_specs{p,4}(1)  : (par_specs{p,4}(end)-par_specs{p,4}(1))/(parsteps-1)  :    par_specs{p,4}(end);

    end
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
    else
        error(['Cannot find requested model in f_modelsettings: ' whichmodels{m}])
    end
end

end

