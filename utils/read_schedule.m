% read a schedule from Wolfgang's discrete planner
% into a [3 x Kpts x Nrobots] waypoint array
%
function [paths, confEllipsoids, obsEllipsoids] = read_schedule(fname)
	json = loadjson(fname);
	N = length(json.agents);
	k = 0;
	for i=1:N
		k = max(k, length(json.agents{i}.path));
	end
	paths = nan(3,k,N);
    confEllipsoids = nan(N,3);
    obsEllipsoids = nan(N,3);
	for i=1:N
        %get path
		p = json.agents{i}.path;
		len = length(p);
		for j=1:len
			paths(1,j,i) = str2num(p{j}.x);
			paths(2,j,i) = str2num(p{j}.y);
			paths(3,j,i) = str2num(p{j}.z);
		end
		for j=(len+1):k
			paths(:,j,i) = paths(:,len,i);
        end
        %get ellipsoid sizes
        ce = json.agents{i}.conflictEllipsoid;
        confEllipsoids(i,1) = str2num(ce.rx);
        confEllipsoids(i,2) = str2num(ce.ry);
        confEllipsoids(i,3) = str2num(ce.rz);
        oe = json.agents{i}.obstacleEllipsoid;
        obsEllipsoids(i,1) = str2num(oe.rx);
        obsEllipsoids(i,2) = str2num(oe.ry);
        obsEllipsoids(i,3) = str2num(oe.rz);
	end
	assert(~any(isnan(paths(:))));
end
