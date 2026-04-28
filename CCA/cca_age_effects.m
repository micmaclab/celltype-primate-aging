% Code inspired by https://github.com/XihanZhang/human-cellular-func-con/blob/main/scripts/08a_PermCCA_GradientVarbyAllCell.m
clear variables
close all


%% ---------------------------
%  Set Up Directories Dynamically
% ---------------------------
script_dir = fileparts(mfilename('fullpath'));
base_dir   = fileparts(fileparts(script_dir)); 

% Define and create output folders
output_main_dir = fullfile(base_dir, 'celltype-primate-aging', 'CCA', 'output');
output_sens_dir = fullfile(base_dir, 'celltype-primate-aging','CCA', 'sensitivity_analysis');
mkdir(output_main_dir);
mkdir(output_sens_dir);

% Add PermCCA code to path (assumes it is inside an 'external' folder in your repo)
addpath(genpath(fullfile(base_dir, 'external', 'PermCCA-master')));


%% ---------------------------
%  Load Shared Data
% ---------------------------

% 1. Age Effects Vector 
age_effects_path  = fullfile(base_dir, '/celltype-primate-aging/Mixed_Effects_Models', 'regionwise_age_effects_MixedLM.csv');
age_effects_table = readtable(age_effects_path);
age_effects_table = sortrows(age_effects_table, 'region');
metric            = age_effects_table.age_t;

% 2. Cell type abundance
cell_data_path = fullfile(base_dir, 'celltype-primate-aging/data/d99_cell_abundance.csv');
Cell_data      = readtable(cell_data_path, 'ReadRowNames', true);
row_names      = str2double(Cell_data.Properties.RowNames);

% Find the index of the row where the name is 70
row_to_remove = row_names == 70;

% Remove that row
Cell_data(row_to_remove, :) = [];
%Cell_data('region', :) = [];

Celldata_mat   = table2array(Cell_data);
cell_type_list = Cell_data.Properties.VariableNames;
cell_num       = size(Celldata_mat, 2);
nP = 5000;  % Number of permutations

%% =========================================================================
%  PART 1 — CCA on Full Data using PermCCA
%% =========================================================================

[p_sim, r_sim, A_sim, B_sim, U_sim, V_sim] = ...
    permcca(metric, Celldata_mat, nP, [], [], [], []);

% Compute loadings and significance
loadings = zeros(cell_num, 2);
for i = 1:cell_num
    [loadings(i, 1), loadings(i, 2)] = corr(Celldata_mat(:, i), V_sim);
end

% Plot
mdl = fitlm(U_sim, V_sim);
figure('Position', [10 10 1300 400]);

% Bar plot of loadings
subplot(1, 2, 1);
bar(loadings(:, 1));
set(gca, 'XTick', 1:cell_num, 'XTickLabels', cell_type_list);
xtickangle(45);
ylabel('Loading');
title('Cell-type Loadings');
for ii = 1:cell_num
    if loadings(ii, 2) < (0.05 / cell_num)
        y = loadings(ii, 1);
        text(ii, y + sign(y)*0.015, '*', 'HorizontalAlignment', 'center');
    end
end

% Scatter plot of canonical correlation
subplot(1, 2, 2);
scatter(U_sim, V_sim); hold on;
plot(mdl);
xlabel('U'); ylabel('V');
title(sprintf('Age Effects r=%.3f, p=%.4f (age_effects)', r_sim, p_sim));

% Save figure
saveas(gcf, fullfile(output_main_dir, 'Cell_age_effects.png'));

% Save main outputs
save(fullfile(output_main_dir, 'permcca_age_effects_Celltypes.mat'), ...
    'p_sim', 'r_sim', 'A_sim', 'B_sim', 'U_sim', 'V_sim', ...
    'loadings', 'cell_type_list');

% Write tables
ColNames = {'loading', 'p_value'};
writetable(array2table(loadings, 'RowNames', cell_type_list, 'VariableNames', ColNames), ...
    fullfile(output_main_dir, 'CCA_loadings_age_effects.csv'), 'WriteRowNames', true);

writetable(array2table([metric, V_sim, U_sim], ...
    'VariableNames', {'age_effects', 'CV_cell', 'CV_age_effects'}), ...
    fullfile(output_main_dir, 'canonical_variables_age_effects.csv'));

writetable(array2table([r_sim, p_sim], ...
    'VariableNames', {'corr_age_effects', 'p_spin_age_effects'}), ...
    fullfile(output_main_dir, 'CCA_significance_age_effects.csv'));

%% =========================================================================
%  PART 2 — Sensitivity Analysis: Gradual Removal of Cell Types
%% =========================================================================

% Define indices
id_L4      = 12;
id_OPC     = 5;
id_ASC     = 1;  % New Index included
id_PV_CHC  = 19;
id_RELN    = 20;

% Define cell type combinations to remove
removal_sets = {
    % Single Sets
    [id_L4],                                    'L4';
    [id_OPC],                                   'OPC';
    [id_ASC],                                   'ASC'; % Added single
    [id_PV_CHC],                                'PV_CHC';
    [id_RELN],                                  'RELN';
    
    % Double Sets
    [id_L4, id_OPC],                            'L4/OPC';
    [id_L4, id_PV_CHC],                         'L4/PV_CHC';
    [id_L4, id_RELN],                           'L4/RELN';
    [id_L4, id_ASC],                            'L4/ASC'; % Added combo
    [id_OPC, id_PV_CHC],                        'OPC/PV_CHC';
    [id_OPC, id_RELN],                          'OPC/RELN';
    [id_OPC, id_ASC],                           'OPC/ASC'; % Added combo
    [id_PV_CHC, id_RELN],                       'PV_CHC/RELN';
    [id_PV_CHC, id_ASC],                        'PV_CHC/ASC'; % Added combo
    [id_RELN, id_ASC],                          'RELN/ASC'; % Added combo
    
    % Triple Sets
    [id_L4, id_OPC, id_PV_CHC],                 'L4/OPC/PV_CHC';
    [id_L4, id_OPC, id_RELN],                   'L4/OPC/RELN';
    [id_L4, id_PV_CHC, id_RELN],                'L4/PV_CHC/RELN';
    [id_OPC, id_PV_CHC, id_RELN],               'OPC/PV_CHC/RELN';
    [id_L4, id_OPC, id_ASC],                    'L4/OPC/ASC'; % Added combo
    
    % Quad & Penta Sets
    [id_L4, id_OPC, id_PV_CHC, id_RELN],        'L4/OPC/PV_CHC/RELN';
    [id_L4, id_OPC, id_PV_CHC, id_RELN, id_ASC], 'L4/OPC/PV_CHC/RELN/ASC'; % Added full set
};

cca_r_p = zeros(length(removal_sets), 2);
RowNames = cell(length(removal_sets), 1);
V_save   = [];

for i = 1:length(removal_sets)
    indices_to_remove = removal_sets{i, 1};
    label             = removal_sets{i, 2};

    Celldata_mat_reduced = Celldata_mat;
    Celldata_mat_reduced(:, indices_to_remove) = [];

    [p_val, r_val, ~, ~, ~, V] = permcca(metric, Celldata_mat_reduced, nP, [], [], [], []);

    cca_r_p(i, :) = [r_val, p_val];
    RowNames{i}   = label;
    V_save        = [V_save, V];
end

cca_table = array2table(cca_r_p, ...
    'RowNames', RowNames, ...
    'VariableNames', {'r', 'p'});
writetable(cca_table, ...
    fullfile(output_sens_dir, 'age_effects_ccaGradualRemove_r_p.csv'), ...
    'WriteRowNames', true);

% --- Output Table 2: Canonical Variables (Regions) ---
V_table = array2table(V_save, ...
    'VariableNames', RowNames);

region_col = table(age_effects_table.region, 'VariableNames', {'region'});
V_table = [region_col, V_table]; 

% 3. Save the table
% Since 'region' is now a standard data column, we don't need 'WriteRowNames'
writetable(V_table, ...
    fullfile(output_sens_dir, 'sensitivity_canonical_variable_age_effects.csv'));
