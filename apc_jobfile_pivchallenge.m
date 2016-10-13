
function JOB_LIST = apc_jobfile_pivchallenge()


region_width  = 128;
region_height = 128 ;

% Decide where the image repo is
if ismac
    img_repo = '~/Desktop/piv_test_images';
%       img_repo = '/Volumes/duo/piv_test_images/';
elseif islinux
    img_repo = '/home/shannon/b/aether/piv_test_images';
end

% case_name = sprintf('Aall');
% JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
% JobFile.Images.BaseName = sprintf('%s_', case_name);
% JobFile.Images.BaseName = 'A';
% JobFile.Images.Extension = '.tif';

case_name = sprintf('Ball');
JobFile.Images.Directory = fullfile(img_repo, case_name, 'raw');
JobFile.Images.BaseName = sprintf('%s_', case_name);
JobFile.Images.BaseName = 'B';
JobFile.Images.Extension = '.bmp';
JobFile.Images.NumDigits = 3;
JobFile.Images.Start = 1;
JobFile.Images.End = 100;
JobFile.Images.Skip = 1;
JobFile.Images.Trailer_01 = 'a';
JobFile.Images.Trailer_02 = 'b';

JobFile.Parameters.Processing.CorrelationStep = 0;
JobFile.Parameters.Processing.RegionWidth  = region_width;
JobFile.Parameters.Processing.RegionHeight = region_height;
JobFile.Parameters.Processing.WindowFraction = 0.5;
JobFile.Parameters.Processing.Shuffle.Range.X = 0;
JobFile.Parameters.Processing.Shuffle.Range.Y = 0;
JobFile.Parameters.Processing.Shuffle.Step.X = 0;
JobFile.Parameters.Processing.Shuffle.Step.Y = 0;

JobFile.Parameters.Processing.Grid.Spacing.X = 16;
JobFile.Parameters.Processing.Grid.Spacing.Y = 16;
JobFile.Parameters.Processing.Grid.Buffer.X = ...
    round(region_width  /2) * [1, 1];
JobFile.Parameters.Processing.Grid.Buffer.Y = ...
    round(region_height /2) * [1, 1];

% RPC Diameter
JobFile.Parameters.Processing.RpcDiameter = 3;

JobFile.JobOptions.DoCorrelations = true;
JobFile.JobOptions.CalculateDisplacements = true;
JobFile.JobOptions.CalculateError = false;

JOB_LIST(1) = JobFile;



end



