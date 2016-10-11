function host_names_list = get_ecn_host_names(host_names_path)

% Default to the list being in the current directory and
% and having this specific name.
if nargin == 0 || isempty(host_names_path)
    host_names_path = 'ecn_host_names';
end

% Open the file
fid = fopen(host_names_path, 'r');

host_names_list = {};
n_hosts = 0;

tline = fgetl(fid);
while ischar(tline)
    if tline(1) ~= '#'
        n_hosts = n_hosts + 1;
        host_names_list{n_hosts} = tline;
    end
    tline = fgetl(fid);
end
fclose(fid);

% % Save output as unique names to prevent replication
host_names_list = (host_names_list)';

end