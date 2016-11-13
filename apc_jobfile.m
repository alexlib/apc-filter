
function JOB_LIST = apc_jobfile()

% Region sizes
region_width  = 128;
region_height = 128 ;

% Directory to the images
JobFile.Images.Directory = '/Users/matthewgiarra/Documents/School/VT/Research/Aether/piv_test_images/piv_challenge_2014_A/raw';

% Image base name
JobFile.Images.BaseName = 'A_';
JobFile.Images.Extension = '.tif';
JobFile.Images.NumDigits = 5;
JobFile.Images.Start = 1;
JobFile.Images.End = 10;
JobFile.Images.Skip = 2;
JobFile.Images.Trailer_01 = 'a';
JobFile.Images.Trailer_02 = 'b';

% Vector (output) paths, names etc
JobFile.Vector.Directory = fullfile(JobFile.Images.Directory, '..', 'vect');
JobFile.Vector.BaseName = JobFile.Images.BaseName;
JobFile.Vector.NumDigits = JobFile.Images.NumDigits;

% Processing parameters (correlation)
JobFile.Parameters.Processing.CorrelationStep = 0;
JobFile.Parameters.Processing.RegionWidth  = region_width;
JobFile.Parameters.Processing.RegionHeight = region_height;
JobFile.Parameters.Processing.WindowFraction = 0.5;


% Grid parameters
JobFile.Parameters.Processing.Grid.Spacing.X = 64;
JobFile.Parameters.Processing.Grid.Spacing.Y = 64;
JobFile.Parameters.Processing.Grid.Buffer.X = round(region_width  /2) * [1, 1];
JobFile.Parameters.Processing.Grid.Buffer.Y = round(region_height /2) * [1, 1];

% RPC diameter
JobFile.Parameters.Processing.RpcDiameter = 3;

% Convergence fraction
JobFile.Parameters.Processing.ConvergenceFraction = 0.01;

% Job Options
JobFile.JobOptions.LoadFilters = false;
JobFile.JobOptions.MakePlots = false;
JobFile.JobOptions.DoCorrelations = false;
JobFile.JobOptions.CalculateDisplacements = false;
JobFile.JobOptions.DoApcFilter = true;

JobFile.Processing.APC.Scheme = 'static';
JobFile.Parameters.Processing.APC.Shuffle.Range.X = 0;
JobFile.Parameters.Processing.APC.Shuffle.Range.Y = 0;
JobFile.Parameters.Processing.APC.Shuffle.Range.Z = 0;
JobFile.Parameters.Processing.APC.Shuffle.Step.X = 0;
JobFile.Parameters.Processing.APC.Shuffle.Step.Y = 0;
JobFile.Parameters.Processing.APC.Shuffle.Step.Z = 0;
JobFile.Parameters.Processing.APC.DiameterUpperBound = 3;

% Append to list
JOB_LIST(1) = JobFile;


end



