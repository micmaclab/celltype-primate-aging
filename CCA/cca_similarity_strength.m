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

% 1. Load and Sort Similarity Strength
sim_path = fullfile(base_dir,'celltype-primate-aging', 'MIND_Network', 'total_similarity_strength.csv');
similarity_strength_table = readtable(sim_path, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');

% Convert RowNames to double for numeric sorting
row_nums = str2double(similarity_strength_table.Properties.RowNames);
[~, sort_idx] = sort(row_nums);
similarity_strength_table = similarity_strength_table(sort_idx, :);

% Extract the metric
metric = similarity_strength_table.("total_similarity_strength");

% Permuted Metrics
%nulls = readmatrix(fullfile(base_dir, 'null_maps_similarity_strength.csv'));
%disp(size(nulls)); 

% Cell type abundance
cell_data_path = fullfile(base_dir,'celltype-primate-aging','data', 'd99_cell_abundance.csv');
Cell_data      = readtable(cell_data_path, 'ReadRowNames', true);
row_names      = str2double(Cell_data.Properties.RowNames);

% Find the index of the row where the name is 70
row_to_remove = row_names == 70;

% Remove that row
Cell_data(row_to_remove, :) = [];

Celldata_mat   = table2array(Cell_data);
cell_type_list = Cell_data.Properties.VariableNames;
cell_num       = size(Celldata_mat, 2);

nP = 5000;  % Number of permutations

%% =========================================================================
%  PART 1 — CCA on Full Data
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
title(sprintf('Similarity Strength r=%.3f, p=%.4f (similarity_strength)', r_sim, p_sim));

% Save figure
saveas(gcf, fullfile(output_main_dir, 'Cell_similarity_strength.png'));

% Save main outputs
save(fullfile(output_main_dir, 'permcca_similarity_strength_Celltypes.mat'), ...
    'p_sim', 'r_sim', 'A_sim', 'B_sim', 'U_sim', 'V_sim', ...
    'loadings', 'cell_type_list');

% Write tables
ColNames = {'loading', 'p_value'};
writetable(array2table(loadings, 'RowNames', cell_type_list, 'VariableNames', ColNames), ...
    fullfile(output_main_dir, 'CCA_loadings_similarity_strength.csv'), 'WriteRowNames', true);

writetable(array2table([metric, V_sim, U_sim], ...
    'VariableNames', {'similarity_strength', 'CV_cell', 'CV_similarity_strength'}), ...
    fullfile(output_main_dir, 'canonical_variables_similarity_strength.csv'));

writetable(array2table([r_sim, p_sim], ...
    'VariableNames', {'corr_similarity_strength', 'p_spin_similarity_strength'}), ...
    fullfile(output_main_dir, 'CCA_significance_similarity_strength.csv'));

%% =========================================================================
%  PART 2 — Sensitivity Analysis: Gradual Removal of Cell Types
%% =========================================================================

% Define indices
id_OLG      = 4;
id_L4_5     = 13;
id_L5_6     = 15;
id_VIP_RELN = 23;
id_L4_5_6   = 14; % New Index included

% Define cell type combinations to remove
removal_sets = {
    % Single Sets
    [id_OLG],                                  'OLG';
    [id_L4_5],                                 'L4_5';
    [id_L5_6],                                 'L5_6';
    [id_L4_5_6],                               'L4_5_6'; % Added single
    [id_VIP_RELN],                             'VIP_RELN';
    
    % Double Sets
    [id_OLG, id_L4_5],                         'OLG/L4_5';
    [id_OLG, id_L5_6],                         'OLG/L5_6';
    [id_OLG, id_L4_5_6],                       'OLG/L4_5_6'; % Added combo
    [id_OLG, id_VIP_RELN],                     'OLG/VIP_RELN';
    [id_L4_5, id_L5_6],                        'L4_5/L5_6';
    [id_L4_5, id_VIP_RELN],                    'L4_5/VIP_RELN';
    [id_L5_6, id_VIP_RELN],                    'L5_6/VIP_RELN';
    [id_L4_5_6, id_VIP_RELN],                  'L4_5_6/VIP_RELN'; % Added combo
    [id_VIP_RELN, id_L5_6],                    'VIP_RELN/L5_6';
    
    % Triple Sets
    [id_OLG, id_L4_5, id_L5_6],                'OLG/L4_5/L5_6';
    [id_OLG, id_L4_5, id_VIP_RELN],            'OLG/L4_5/VIP_RELN';
    [id_OLG, id_L5_6, id_VIP_RELN],            'OLG/L5_6/VIP_RELN';
    [id_OLG, id_L4_5_6, id_VIP_RELN],          'OLG/L4_5_6/VIP_RELN'; % Added combo
    [id_L4_5, id_L5_6, id_VIP_RELN],           'L4_5/L5_6/VIP_RELN';
    
    % Quad Sets
    [id_OLG, id_L4_5, id_L5_6, id_VIP_RELN],   'OLG/L4_5/L5_6/VIP_RELN';
    [id_OLG, id_L4_5, id_L5_6, id_L4_5_6],   'OLG/L4_5/L5_6/L4_5_6';
    [id_OLG, id_L4_5, id_L5_6, id_L4_5_6, id_VIP_RELN],   'OLG/L4_5/L5_6/L4_5_6/VIP_RELN';
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
region_ids = similarity_strength_table.Properties.RowNames;

% 2. Create the cca_table (Statistics)
cca_table = array2table(cca_r_p, ...
    'RowNames', RowNames, ...
    'VariableNames', {'r', 'p'});

% 3. Create the V_table (Canonical Variables)
V_table = array2table(V_save, ...
    'VariableNames', RowNames);

% 4. Create the region column and prepend it to V_table
% We convert region_ids to a table first, then concatenate
region_col = table(region_ids, 'VariableNames', {'region'});
V_table = [region_col, V_table]; 

% 5. Save results
% Use 'WriteRowNames' for cca_table because the labels are in the RowNames property
writetable(cca_table, ...
    fullfile(output_sens_dir, 'similarity_strength_ccaGradualRemove_r_p.csv'), ...
    'WriteRowNames', true);

% For V_table, 'region' is now a standard column, so 'WriteRowNames' is not needed
writetable(V_table, ...
    fullfile(output_sens_dir, 'sensitivity_canonical_variable_similarity_strength.csv'));
