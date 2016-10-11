function [N, S] = par_hosts(host_names_path)

% Default to not inputting a host name list path
if nargin < 1
    host_names_path = [];
end

% Get the list of host names
host_names_list = (get_ecn_host_names(host_names_path));

% Number of hosts
N = length(host_names_list);

% Get the name of this machine
[~, my_name] = system('hostname');

% Line feed character in ascii
lfchar = char(10);

% truncate the EOL character
my_name = strrep(my_name, lfchar, '');

% Find your position within the list of hosts
S = find(strcmp(my_name, host_names_list));

end