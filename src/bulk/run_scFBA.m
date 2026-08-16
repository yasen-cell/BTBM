%% scFBA pipeline for bulk data (HGNC_ID-ready CSV)
clear; clc;

%% 1. Load bulk transcriptomic data (e.g., TNBC, MTC84, NCI60 ,or ccRCC)
fprintf('Reading input CSV...\n');
T = readtable('bulk_data.csv', 'VariableNamingRule', 'preserve'); 

%% 2. Load HMRcore model
fprintf('Loading HMRcore model...\n');
tmp = load('HMRcore.mat');
HMRcore = tmp.HMRcore;

%% 3. Identify required columns
varNames = T.Properties.VariableNames;
if ~ismember('HGNC_ID', varNames), error('Column HGNC_ID not found.'); end
if ~ismember('Pooled', varNames), error('Column Pooled not found.'); end

sampleMask = startsWith(varNames, 'GSM');
SampleNames = varNames(sampleMask);
if isempty(SampleNames), error('No GSM sample columns detected.'); end

%% 4. Clean HGNC_ID and align to HMRcore gene order
rawHGNC = strtrim(string(T.HGNC_ID));
validIdx = ~(ismissing(rawHGNC) | rawHGNC == "");
T = T(validIdx, :);
rawHGNC = rawHGNC(validIdx);

[rawHGNC_unique, ia] = unique(rawHGNC, 'stable');
T_unique = T(ia, :);

modelGenes = strtrim(string(HMRcore.genes));
[tf, loc] = ismember(modelGenes, rawHGNC_unique);

BC04_filt = T_unique(loc(tf), :);
BC04_filt.HGNC_ID = cellstr(modelGenes(tf));

if height(BC04_filt) == 0, error('No genes matched HMRcore after alignment.'); end

%% 5. Build expression matrices
pooledExpr = BC04_filt.Pooled;
exprMat = BC04_filt{:, SampleNames};

if any(isnan(exprMat(:))) || any(isnan(pooledExpr(:)))
    error('Expression data contains NaN values.');
end

%% 6. Build dataset structure
geneList = strtrim(cellstr(BC04_filt.HGNC_ID));
BC04 = makeSCdataset(pooledExpr, exprMat, SampleNames, geneList, 1e-4);
BC04 = Genes_Sign(BC04);
BC04 = RepairNegFalse(BC04);

%% 7. Set solver to Gurobi
changeCobraSolver('gurobi');
if ~strcmpi(getCobraSolver('LP'), 'gurobi')
    error('Failed to set Gurobi solver.');
end

%% 8. Build integrated population model
[~, Ex_id] = EditBoundaries(HMRcore, 'Ex_');
IdxExRxns = Ex_id.ID;

[~, Coop_id] = EditBoundaries(HMRcore, '_COOP');
[~, Biomass_id] = EditBoundaries(HMRcore, 'biomass_synthesis');
IdxCoopRxn = [Coop_id.ID; Biomass_id.ID];

HMRcore = ScFBAExpSetting(HMRcore, length(BC04.CellType));
BC04 = single2IntPopModel(BC04, HMRcore, IdxExRxns, IdxCoopRxn, 's');

%% 9. Optimize integrated model
fluxIntBC04 = optimizeCbModel(BC04.modelFVAInt);
if isempty(fluxIntBC04.x), error('Optimization failed: no flux solution.'); end

%% 10. Split fluxes and export
[FluxPopBC04, RxnPop] = splitScFluxes(BC04.modelFVAInt, fluxIntBC04, length(BC04.CellType));

FluxTable = array2table(FluxPopBC04, 'VariableNames', SampleNames);
FluxTable = addvars(FluxTable, RxnPop, 'Before', 1, 'NewVariableNames', 'Reaction');

writetable(FluxTable, 'FluxTable_scFBA.csv');
fprintf('Done. Results saved to FluxTable_scFBA.csv\n');