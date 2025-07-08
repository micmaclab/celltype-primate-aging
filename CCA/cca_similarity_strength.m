% Code inspired by https://github.com/XihanZhang/human-cellular-func-con/blob/main/scripts/08a_PermCCA_GradientVarbyAllCell.m
clear variables
close all

%% ---------------------------
%  Set Up Directories Dynamically
% ---------------------------

this_file  = matlab.desktop.editor.getActiveFilename;
script_dir = fileparts(this_file);
base_dir   = fileparts(script_dir);  % Assumes script is in home/CCA/

% Define subdirectories
data_dir         = fullfile(base_dir, 'CCA', 'input');
output_main_dir  = fullfile(base_dir, 'CCA', 'output');
output_sens_dir  = fullfile(base_dir, 'CCA', 'sensitivity_analysis');

% Create output folders if they don't exist
if ~exist(output_main_dir, 'dir'), mkdir(output_main_dir); end
if ~exist(output_sens_dir, 'dir'), mkdir(output_sens_dir); end

% Add PermCCA code to path - change to correct path as needed 
addpath('/Users/melinatsotras/Downloads/PermCCA-master');

%% ---------------------------
%  Load Shared Data
% ---------------------------

% 1. Total Similarity Strength
similarity_strength_table = readtable(fullfile(base_dir,'MIND_Network', 'total_similarity_strength.csv'), 'ReadRowNames', true);

% Get the row names (cell array of char)
row_names = similarity_strength_table.Properties.RowNames;

% Convert them to numeric
row_nums = str2double(row_names);

% Sort numeric row names and get sorting indices
[sorted_nums, sort_idx] = sort(row_nums);

% Reorder the table rows according to sort order
similarity_strength_table = similarity_strength_table(sort_idx, :);

% Update row names to sorted (optional but cleaner)
similarity_strength_table.Properties.RowNames = string(sorted_nums);


% Extract the 't_value' column as a numeric vector
metric = similarity_strength_table.total_similarity_strength;

% 2. Permuted Metrics
nulls = readmatrix(fullfile(data_dir, 'hungarian_5k_nulls_similarity_strength.csv'));
disp(size(nulls));  % Confirm dimensions

% 3. Cell type abundance
Cell_data      = readtable(fullfile(base_dir, 'data', 'd99_cell_abundance.csv'), 'ReadRowNames', true);
row_names = str2double(Cell_data.Properties.RowNames);

% Find the index of the row where the name is 70
row_to_remove = row_names == 70;

% Remove that row
Cell_data(row_to_remove, :) = [];

Cell_data{:, :} = Cell_data{:, :} / 100000;
Celldata_mat   = table2array(Cell_data);
cell_type_list = Cell_data.Properties.VariableNames;
cell_num       = size(Celldata_mat, 2);

nP = 5000;  % Number of permutations

%% =========================================================================
%  PART 1 — CCA on Full Data
%% =========================================================================

[p_sim, r_sim, A_sim, B_sim, U_sim, V_sim] = ...
    permcca(metric, Celldata_mat, nP, [], [], [], [], nulls);

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
id_OLG  = 4;
id_L4_5 = 13;
id_L5_6 = 15;
id_RELN = 20;

% Define cell type combinations to remove
removal_sets = {
    [id_OLG, id_L4_5, id_L5_6],            'OLG/L4_5/L5_6';
    [id_L4_5],                             'L4_5';
    [id_OLG],                              'OLG';
    [id_RELN],                             'RELN';
    [id_L5_6],                             'L5_6';
    [id_L4_5, id_RELN],                    'L4_5/RELN';
    [id_OLG, id_L5_6],                     'OLG/L5_6';
    [id_L4_5, id_L5_6],                    'L4_5/L5_6';
    [id_OLG, id_L4_5],                     'OLG/L4_5';
    [id_RELN, id_L5_6],                    'RELN/L5_6';
    [id_OLG, id_L4_5, id_L5_6, id_RELN],   'OLG/L4_5/RELN/L5_6';
    [id_OLG, id_RELN],                     'OLG/RELN';
    [id_L5_6, id_RELN],                    'L5_6/RELN';
    [id_OLG, id_L4_5, id_RELN],            'OLG/L4_5/RELN';
    [id_OLG, id_L5_6, id_RELN],            'OLG/L5_6/RELN';
    [id_L4_5, id_L5_6, id_RELN],           'L4_5/L5_6/RELN';
};

cca_r_p = zeros(length(removal_sets), 2);
RowNames = cell(length(removal_sets), 1);
V_save   = [];

for i = 1:length(removal_sets)
    indices_to_remove = removal_sets{i, 1};
    label             = removal_sets{i, 2};

    Celldata_mat_reduced = Celldata_mat;
    Celldata_mat_reduced(:, indices_to_remove) = [];

    [p_val, r_val, ~, ~, ~, V] = permcca(metric, Celldata_mat_reduced, nP, [], [], [], [], nulls);

    cca_r_p(i, :) = [r_val, p_val];
    RowNames{i}   = label;
    V_save        = [V_save, V];
end

% Output tables
cca_table = array2table(cca_r_p, ...
    'RowNames', RowNames, ...
    'VariableNames', {'r', 'p'});

V_table = array2table(V_save, ...
    'VariableNames', matlab.lang.makeValidName(RowNames));

% Save results
writetable(cca_table, ...
    fullfile(output_sens_dir, 'similarity_strength_ccaGradualRemove_r_p.csv'), ...
    'WriteRowNames', true);

writetable(V_table, ...
    fullfile(output_sens_dir, 'sensitivity_canonical_variable_similarity_strength.csv'));
