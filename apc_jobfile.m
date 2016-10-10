
% Diffusion standard deviation
diffusion_std = 3;

region_width  = 128;
region_height = 128 ;

img_repo = '/Users/matthewgiarra/Desktop/piv_test_images';
case_name = sprintf('poiseuille_diffusion_%0.2f', diffusion_std);


solution_file_name = sprintf('poiseuille_diffusion_%0.2f.mat', diffusion_std);


JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
JobFile.Solution.Path = fullfile(img_repo, case_name, 'jobfiles', solution_file_name);

JobFile.Images.BaseName = sprintf('%s_', case_name);
JobFile.Images.Extension = '.tiff';
JobFile.Images.NumDigits = 6;
JobFile.Images.Start = 1;
JobFile.Images.End = 50;
JobFile.Images.Skip = 2;
JobFile.Images.Trailer_01 = '';
JobFile.Images.Trailer_02 = '';

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