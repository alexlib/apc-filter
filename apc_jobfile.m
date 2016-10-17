
function JOB_LIST = apc_jobfile()

region_width  = 128;
region_height = 128 ;

% Decide where the image repo is
img_repo = get_image_repo();

% Diffusion standard deviation
diffusion_std = 1;

case_name = sprintf('poiseuille_diffusion_%0.2f', diffusion_std);


solution_file_name = sprintf('poiseuille_diffusion_%0.2f.mat', diffusion_std);


JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
JobFile.Solution.Path = fullfile(img_repo, case_name, 'jobfiles', solution_file_name);

JobFile.Images.BaseName = sprintf('%s_', case_name);
JobFile.Images.Extension = '.tiff';
JobFile.Images.NumDigits = 6;
JobFile.Images.Start = 1;
JobFile.Images.End = 2000;
JobFile.Images.Skip = 2;
JobFile.Images.Trailer_01 = '';
JobFile.Images.Trailer_02 = '';
JobFile.Images.Height = 2048;
JobFile.Images.Width = 2048;

JobFile.Parameters.Processing.CorrelationStep = 1;
JobFile.Parameters.Processing.RegionWidth  = region_width;
JobFile.Parameters.Processing.RegionHeight = region_height;
JobFile.Parameters.Processing.WindowFraction = 0.5;
JobFile.Parameters.Processing.Shuffle.Range.X = 0;
JobFile.Parameters.Processing.Shuffle.Range.Y = 0;
JobFile.Parameters.Processing.Shuffle.Step.X = 0;
JobFile.Parameters.Processing.Shuffle.Step.Y = 0;

JobFile.Parameters.Processing.Grid.Spacing.X = 64;
JobFile.Parameters.Processing.Grid.Spacing.Y = 64;
JobFile.Parameters.Processing.Grid.Buffer.X = round(region_width  /2) * [1, 1];
JobFile.Parameters.Processing.Grid.Buffer.Y = round(region_height /2) * [1, 1];
% JobFile.Parameters.Processing.Grid.Buffer.X = round(JobFile.Images.Width / 2) * [1, 1];
% JobFile.Parameters.Processing.Grid.Buffer.Y = round(JobFile.Images.Height / 2) * [1, 1];


JobFile.Parameters.Processing.RpcDiameter = 3;

JobFile.JobOptions.LoadFilters = true;

JobFile.JobOptions.DoCorrelations = false;
JobFile.JobOptions.CalculateDisplacements = false;
JobFile.JobOptions.CalculateError = true;


JOB_LIST(1) = JobFile;

% % New Std Dev
diffusion_std = 3;
case_name = sprintf('poiseuille_diffusion_%0.2f', diffusion_std);
solution_file_name = sprintf('poiseuille_diffusion_%0.2f.mat', diffusion_std);
JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
JobFile.Solution.Path = fullfile(img_repo, case_name, 'jobfiles', solution_file_name);
JobFile.Images.BaseName = sprintf('%s_', case_name);

JOB_LIST(end + 1) = JobFile;
%  
% % % New Std Dev
diffusion_std = 5;
case_name = sprintf('poiseuille_diffusion_%0.2f', diffusion_std);
solution_file_name = sprintf('poiseuille_diffusion_%0.2f.mat', diffusion_std);
JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
JobFile.Solution.Path = fullfile(img_repo, case_name, 'jobfiles', solution_file_name);
JobFile.Images.BaseName = sprintf('%s_', case_name);
% % 
JOB_LIST(end + 1) = JobFile;


end



