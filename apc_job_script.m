function apc_job_script()

% Add paths.
addpaths('..');

% Get your number from the list of hosts
[N, S] = par_hosts();

% Read the job file;
JobList = apc_jobfile;

% NUmber of jobs
num_jobs = length(JobList);

start_pair = 1;
end_pair = 1000;

for n = 1 : num_jobs
    
    JobFile = JobList(n);
    
    c_step = JobFile.Parameters.Processing.CorrelationStep;
    
    % Read the image step
    img_skip = JobFile.Images.Skip;
    
    % Sets vector
    pairs_vect = round(linspace(start_pair, end_pair + 1, N+1));

    % Start and end sets
    start_pair_current = pairs_vect(S);
    end_pair_current = pairs_vect(S + 1) - 1;

    start_image_current = (start_pair_current - 1) * img_skip + 1;
    end_image_current = (end_pair_current - 1) * img_skip + 1 + c_step;

    % Update the jobfile
    JobList(n).Images.Start = start_image_current;
    JobList(n).Images.End = end_image_current;
    
end

% Run everything
apc_error_analysis_mc_full_images_saveplanes(JobList);

end