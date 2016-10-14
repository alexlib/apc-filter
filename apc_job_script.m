function apc_job_script(do_gen, do_corr, do_fit)

% Add paths.
addpaths('..');

% Get your number from the list of hosts
[N, S] = par_hosts();

% Read the job file;
JobList = apc_jobfile;

% Number of jobs
num_jobs = length(JobList);

% dx_rand_list = [1, 3, 5];
dx_rand_list = 3;
num_pairs = 1000;
image_size = [2048, 2048];

% Decide where the image repo is
img_repo = get_image_repo();

% if nargin < 3
%     job_nums = 1 : num_jobs;
% end

% % Which pairs to do broh
% start_pair = 1;
% end_pair = 1000;

% 
% 
% % Process the data broh.
% for n = 1 : num_jobs
%     
%     % Pick out the job for this one, broh. 
%     JobFile = JobList(n);
%     
%     % This is the correlation step, broh.
%     % Make sure to get this one right, broh.
%     c_step = JobFile.Parameters.Processing.CorrelationStep;
%     
%     % Read the image step
%     img_skip = JobFile.Images.Skip;
%     
%     % Sets vector
%     pairs_vect = round(linspace(start_pair, end_pair + 1, N+1));
% 
%     % Start and end sets
%     start_pair_current = pairs_vect(S);
%     end_pair_current = pairs_vect(S + 1) - 1;
% 
%     start_image_current = (start_pair_current - 1) * img_skip + 1;
%     end_image_current = (end_pair_current - 1) * img_skip + 1 + c_step;
% 
%     % Update the jobfile
%     JobList(n).Images.Start = start_image_current;
%     JobList(n).Images.End = end_image_current;
%     
% end

% Only run if the job is within scope
if S <= num_jobs
    
    if do_gen
        make_poiseuille_images(dx_rand_list(S), ...
            num_pairs, image_size, img_repo);
    end

    % Run everything
    if do_corr
        apc_error_analysis_mc_full_images_saveplanes(JobList(S));
    end
% 
%     % Displacement estimate.
    if do_fit
        apc_error_analysis_mc_full_images_loadplanes(JobList(S));
    end

end



end