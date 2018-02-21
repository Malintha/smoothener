% read a schedule from Wolfgang's discrete planner
% into a [3 x Kpts x Nrobots] waypoint array
%
function [paths,names,types] = read_schedule(fname, read_agents)
	json = loadjson(fname);
 
    if ~exist('read_agents','var')
        read_agents = 1:length(json.agents);
    end
    
    N = length(read_agents);
 
	k = 0;
	for i=read_agents
		k = max(k, length(json.agents{i}.path));
	end
	paths = nan(3,k,N);
    names = cell(N,1);
    types = cell(N,1);
    iprime = 1;
	for i=read_agents
        %get path
		p = json.agents{i}.path;
		len = length(p);
        for j=1:len
            paths(1,j,iprime) = str2num(p{j}.x);
            paths(2,j,iprime) = str2num(p{j}.y);
            paths(3,j,iprime) = str2num(p{j}.z);
        end
        
        for j=(len+1):k
            paths(:,j,iprime) = paths(:,len,iprime);
        end
        %get name
        names{iprime} = json.agents{i}.name;
        %get type
        types{iprime} = json.agents{i}.type;
        
        iprime = iprime + 1;
	end
	assert(~any(isnan(paths(:))));
end
