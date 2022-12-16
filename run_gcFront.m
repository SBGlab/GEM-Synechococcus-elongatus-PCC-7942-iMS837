initCobraToolbox(false)
changeCobraSolver('gurobi')

%SETTING PROJECT PATH:
PROJECT_PATH='/home/alvaro/github/ConsortiumEngineeringTool/jupyter_notebook/'

%READING 'production_target_tf.csv' TO LOAD TARGET METABOLITE & ORGANISM
targetTable=readtable(strcat(PROJECT_PATH,'production_target_tf.csv'))
experiment=char(table2cell(targetTable(:,1)))
%The reaction/metabolite in the model to be growth-coupled.
%If a metabolite name is supplied, then the exchange reaction for this metabolite will be selected
target=char(table2cell(targetTable(:,2)))
framework=char(table2cell(targetTable(:,3)))
biomass_limit=cell2mat(table2cell(targetTable(:,4)))

%SETTING INTERNAL PATHS:
EXPERIMENT_PATH=strcat(PROJECT_PATH, experiment,'_ko_results/')
FRAMEWORK_PATH=strcat(EXPERIMENT_PATH,framework,'/')

%READING A FILE WITH THE NEW MODEL, CONFIGURED WITH Cameo
MODEL_PATH=strcat(EXPERIMENT_PATH,'configured_model.mat')
model=readCbModel(MODEL_PATH)

%READING '{experiment}_optknock_analysis.csv' IMPORTS THE TOP CANDIDATE GENES
%CHOSEN BY OptKnock ALGORITHM EXECUTED BY Cameo
selectedRxnsTable=readtable(strcat(FRAMEWORK_PATH,experiment,'_optknock_gene_analysis.csv'))
selectedRxns=table2cell(selectedRxnsTable(:,1))
modelRxns = model.rxns
noKOsRxns = setdiff(modelRxns, selectedRxns)

%NOW WE NEED TO READ THE FILE 'workflow_params.csv', WHERE ALL
maxKO=30
minGrowth=biomass_limit
minProduction=10^-3

%EVENTUALLY RUN gcFront ALGORITHM

% set model options so gurobi is used to solve LP problems, and so knockout
% combinations with more than 30 knockouts will not be considered
options=struct;
options.solver='gurobi';           % Character array of the name of the LP solver used for solving FBA problems. Default  =  currently set LP solver. If no solver set, will use whatever solver the function initCobraToolbox sets�
options.maxknockouts=maxKO;        % Double specifying the maximum number of knockouts that a design may contain. Default  =  inf (i.e. no limit on the number of deletions allowed)

% Other options that are not being set here:

%options.biomassrxn='';          % Character array specifying the reaction that represents cell growth. Default  =  current objective reaction in model
options.tol=10^-8;               % Double specifying tolerance to mathematical errors. Flux differences smaller than this will be ignored. Default  =  10^-8
options.tiltval=10^-4;           % Double specifying the coefficient used for objective value tilting. Default  =  10^-4
options.shiftval=10^-5;          % Double specifying the size of the change in growth rate used to calculate shadow price of uncoupled designs. Default  =  10^-5
options.skipreduction=true;     % Logical that determines if algorithm should skip deletion of inactive reactions and pooling of unbranched pathways/redundant genes. Default  =  false (i.e. carry out reduction)
options.mingrowth=minGrowth;     % Double specifying the minimum growth threshold- do not consider any deletion that lower growth below this value. Default  =  10^-3
options.minprod=minProduction;   % Double specifying the minimum product threshold- do not consider any deletion that lowers maxmimum product synthesis below this value. Default  =  10^-3
options.removeredundancy=true;   % Logical that detemines if any redundant deletions should be removed from the designs that the GA identifies. Default  =  true
options.newredundantremoval=true;% Logical specifying whether to use a new, faster method for removing redundant KOs, or to use the original version that was used to obtain the data for the gcFront paper.
options.maxreductionsize=maxKO;  % Double specifying the maximum size of design that should be fully explored for removal of redundant KOs
options.saveresults=true;        % Logical that determines if results and algorithm parameters shoudl be saved. Default  =  true
options.deletegenes=false;       % Logical that determines if algorithm should test gene knockouts or reaction knockouts. Default  =  false (i.e. algorithm will knock out reactions)
options.ignorelistrxns=noKOsRxns;       % cell array of reactions that user does not want to be knocked knocked out (e.g. reactions that are associated with experimentally inaccessible genes). Default  =  {} (i.e. no reactions)
options.ignorelistgenes={};      % cell array of genes that should not be knocked out (e.g. essential genes). If reaction knockouts are being tested, reactions that can't be knocked out without knocking out these genes will be ignored. Default  =  {} (i.e. no genes).
options.dontkoess=true;          % Logical to determine whether reactions that can only be knocked out if an essential gene is knocked out are ignored. Default  =  true (i.e do not consider these reactions)
options.onlykogeneassoc=true;    % Logical to determine if reactions that are not gene associated should be ignored. Default  =  true (i.e. only consider deletion of gene associated reactions)
options.mutationrate=1;          % Double controlling the mean number of mutations that the GA mutation function introduces. Default  = 1
options.popsize=1000;             % Double controlling the number of individuals in the GA population. Default  = 200
options.genlimit=100000;         % Double specifying the maximum number of generations the GA will run for. Default  = 10000
options.timelimit=60*60*24;      % Double specifying the maximum length of time (in seconds) that the GA will run for. Default  =  86400 (i.e. 1 day)
options.fitnesslimit=inf;        % Double. The algorithm will terminate if a design with a product synthesis rate higher than this is discovered. Default = inf (i.e. this condition will not stop the algorithm)
options.spreadchangelimit=10^-4; % Double specifying how low the change in spread must be before the algorithm terminates. For further information see MATLAB's 'gamultiobj Algorithm' documentation. Default  =  10^-4�
options.stallgenlimit=[];        % Double specifying the number of generations that change in spread must be below spreadchangelimit before the algorithm terminates. For further information see MATLAB's 'gamultiobj Algorithm' documentation. Default  =  number of generations (i.e. the algorithm will never terminate due to spread being too low)
options.plotinterval=1;          % Double specifying the number of generations that must pass before the algorithm updates the plot of the current designs. Default  =  1 (i.e. updates every generation)
[designTable,algParams,reducedModel]=gcFront(model,target,options);

writetable(designTable, strcat(EXPERIMENT_PATH,'gcfront_strategies.csv'))
