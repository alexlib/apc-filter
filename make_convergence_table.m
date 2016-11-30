
% List of diffusions to use
diffusion_std_list = [1.5, 3.0, 4.5];

% Parameters for calculations 
% This is the number of elements over which the
% moving mean error is calculated
kernel_length = 25;

% This is the convergence criterion
conv_criterion = 0.05;


% Number of diffusions
num_diffusions = length(diffusion_std_list);

% Allocate arrays for results
conv_err_mean_scc = zeros(3, num_diffusions);
conv_err_mean_rpc = zeros(3, num_diffusions);
conv_err_mean_apc = zeros(3, num_diffusions);
conv_err_std_scc = zeros(3, num_diffusions);
conv_err_std_rpc = zeros(3, num_diffusions);
conv_err_std_apc = zeros(3, num_diffusions);

conv_iter_scc = zeros(num_diffusions, 1);
conv_iter_rpc = zeros(num_diffusions, 1);
conv_iter_apc = zeros(num_diffusions, 1);

for k = 1 : num_diffusions
    
    % Extract the diffusion std deviation
    diffusion_std = diffusion_std_list(k);
   
    % Form the path to the results
    results_path = sprintf('/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/poiseuille_diffusion_%0.2f/error/poiseuille_error_h128_w128_diff_std_%0.2f_000001_002000.mat', diffusion_std, diffusion_std);
    
    % Load the results
    load(results_path);
    
    % Read the particle diameter
    dp = particle_diameter;

    % Calculate the particle diameter ratio
    dp_ratio = diffusion_std / (dp / 4);
    
    % Calculate SCC residuals
    [conv_iter_scc(k), diff_x_scc_mean_smoothed, diff_y_scc_mean_smoothed, image_nums] = ...
        calculate_vector_convergence(tx_err_scc_spatial_mat,...
        ty_err_scc_spatial_mat, conv_criterion, kernel_length);

    % Calculate RPC residuals
    [conv_iter_rpc(k), diff_x_rpc_mean_smoothed, diff_y_rpc_mean_smoothed] = ...
        calculate_vector_convergence(tx_err_rpc_spatial_mat,...
        ty_err_rpc_spatial_mat, conv_criterion, kernel_length);

    % Calculate APC residuals
    [conv_iter_apc(k), diff_x_apc_mean_smoothed, diff_y_apc_mean_smoothed] = ...
        calculate_vector_convergence(tx_err_apc_mat,...
        ty_err_apc_mat, conv_criterion, kernel_length);
    
    % Convergence iterations
    conv_iterations = [conv_iter_apc(k); conv_iter_rpc(k); conv_iter_scc(k)];
   
    % Populate the error arrays.
    conv_err_mean_scc(:, k) = mean_err_scc_spatial(conv_iterations);
    conv_err_mean_rpc(:, k) = mean_err_rpc_spatial(conv_iterations);
    conv_err_mean_apc(:, k) = mean_err_apc(conv_iterations);
    
     % Populate the std dev arrays.
    conv_err_std_scc(:, k) = std_err_scc_spatial(conv_iterations);
    conv_err_std_rpc(:, k) = std_err_rpc_spatial(conv_iterations);
    conv_err_std_apc(:, k) = std_err_apc(conv_iterations);
    
    
end

% Results array
results_array_mean = [conv_err_mean_apc(:)'; conv_err_mean_rpc(:)'; conv_err_mean_scc(:)'];
results_array_std = [conv_err_std_apc(:)'; conv_err_std_rpc(:)'; conv_err_std_scc(:)'];
% 
% % Input data for the latex table
input.data = results_array_mean;

% Set column labels (use empty string for no label):
input.tableColLabels = {'APC','RPC','SCC','APC','RPC','SCC','APC','RPC','SCC'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'APC','RPC','SCC'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.2f',9}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 1;

% Uses booktabs basic formating rules ('1' = using booktabs, '0' = not using booktabs). 
% Note that this option requires the booktabs package being available in your LaTex. 
% Also, setting the booktabs option to '1' overwrites input.tableBorders if it exists.
% input.booktabs = 0;


% LaTex table caption:
input.tableCaption = 'MyTableCaption';

% LaTex table label:
input.tableLabel = 'MyTableLabel';

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% call latexTable:
latex_mean = latexTable(input);

clear 'input'
% % Input data for the latex table
input.data = results_array_std;

% Set column labels (use empty string for no label):
input.tableColLabels = {'APC','RPC','SCC','APC','RPC','SCC','APC','RPC','SCC'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'APC','RPC','SCC'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.2f',9}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 1;

% Uses booktabs basic formating rules ('1' = using booktabs, '0' = not using booktabs). 
% Note that this option requires the booktabs package being available in your LaTex. 
% Also, setting the booktabs option to '1' overwrites input.tableBorders if it exists.
% input.booktabs = 0;


% LaTex table caption:
input.tableCaption = 'MyTableCaption';

% LaTex table label:
input.tableLabel = 'MyTableLabel';

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% call latexTable:
latex_std = latexTable(input);




