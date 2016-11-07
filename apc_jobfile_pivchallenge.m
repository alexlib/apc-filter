
function JOB_LIST = apc_jobfile_pivchallenge()


region_width  = 128;
region_height = 128 ;

img_repo = get_image_repo();

% case_name = sprintf('Aall');
% JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
% JobFile.Images.BaseName = sprintf('%s_', case_name);
% JobFile.Images.BaseName = 'A';
% JobFile.Images.Extension = '.tif';

case_name = sprintf('piv_challenge_2014_A');
JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
JobFile.Images.BaseName = sprintf('%s_', case_name);
JobFile.Images.BaseName = 'A_';
JobFile.Images.Extension = '.tif';
JobFile.Images.NumDigits = 5;
JobFile.Images.Start = 1;
JobFile.Images.End = 100;
JobFile.Images.Skip = 1;
JobFile.Images.Trailer_01 = '_a';
JobFile.Images.Trailer_02 = '_b';

JobFile.Parameters.Processing.CorrelationStep = 0;
JobFile.Parameters.Processing.RegionWidth  = region_width;
JobFile.Parameters.Processing.RegionHeight = region_height;
JobFile.Parameters.Processing.WindowFraction = 0.5;
JobFile.Parameters.Processing.Shuffle.Range.X = 0;
JobFile.Parameters.Processing.Shuffle.Range.Y = 0;
JobFile.Parameters.Processing.Shuffle.Step.X = 0;
JobFile.Parameters.Processing.Shuffle.Step.Y = 0;

JobFile.Parameters.Processing.Grid.Spacing.X = 64;
JobFile.Parameters.Processing.Grid.Spacing.Y = 64;
JobFile.Parameters.Processing.Grid.Buffer.X = ...
    round(region_width  /2) * [1, 1];
JobFile.Parameters.Processing.Grid.Buffer.Y = ...
    round(region_height /2) * [1, 1];

% Convergence fraction
JobFile.Parameters.Processing.ConvergenceFraction = 0.01;

% Load filters?
JobFile.JobOptions.LoadFilters = false;

% Make plots?
JobFile.JobOptions.MakePlots = false;

% RPC Diameter
JobFile.Parameters.Processing.RpcDiameter = 3;

JobFile.JobOptions.DoCorrelations = false;
JobFile.JobOptions.CalculateDisplacements = true;
JobFile.JobOptions.CalculateError = false;

JOB_LIST(1) = JobFile;



end



