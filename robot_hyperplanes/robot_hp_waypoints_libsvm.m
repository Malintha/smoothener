function [ A, b ] = robot_hp_waypoints_libsvm( paths, ellipsoids )
%ROBOT_HP_WAYPOINTS Calculates robot seperating hyperplanes
%INPUT:
%   paths: [#dimensions, #waypoints, #robots]
%          paths(:,:,n) are the waypoints for robot n
%   ellipsoids: [#robots, #dimensions]
%          ellipsoids(n,:) the axis-aligned ellipsoid axes size for
%          robot n
%OUTPUT:
%   A: [DIM x NROBOTS x NROBOTS x (NPTS - 1)] array of 
%      hyperplane normal vectors for each robot-robot interaction
%      at each segment.
%   b: distance from origin for hyperplanes. i.e. a(:,...)^T x <= b
%NOTES:
%   currently only works for 3d paths / ellipsoids

Nrob = size(paths,3); % #robots
Nsteps = size(paths,2) -1; % #time steps

%generate untransformed-ellipsoid vertices for all vertices
ntheta = 32;
nphi = 32;

%unit sphere
sphereVerts = ellipsoid_verts(1,1,1,ntheta,nphi);

%generate rest of ellipsoids by scaling sphere
ellVerts = zeros(size(sphereVerts,1),3,Nrob);
for r = 1:Nrob
    ellVerts(:,:,r) = [sphereVerts(:,1).*ellipsoids(r,1),...
                       sphereVerts(:,2).*ellipsoids(r,2),...
                       sphereVerts(:,3).*ellipsoids(r,3)];
end

 %TODO ABOVE: uniform meshes for ellipsoids. 
 %  http://persson.berkeley.edu/distmesh/

NellVerts = size(ellVerts,1); %number of vertices per ellipsoid

%translate ellipsoids to create 'swept hull' vertices for each timestep
hullVerts = zeros(2*NellVerts,size(ellVerts,2),Nrob, Nsteps);
for r = 1:Nrob
    for step = 1:Nsteps
        t1 = repmat(paths(:,step,r)',NellVerts,1);
        t2 = repmat(paths(:,step+1,r)',NellVerts,1);
        hullVerts(1:NellVerts,:,r,step) = ellVerts(:,:,r) + t1;
        hullVerts((NellVerts+1):end,:,r,step) = ellVerts(:,:,r) + t2;      
    end
end

%labels for svm
%   1 for robot 1, -1 for robot 2
labels = [ones(2*NellVerts,1);-1*ones(2*NellVerts,1)];

%initialize output structures
A = nan(3,Nrob,Nrob,Nsteps);
b = nan(Nrob,Nrob,Nsteps);

%parfor each timestep
parfor step = 1:Nsteps
    %renames for parfor
    stepA = nan(3,Nrob,Nrob);
    stepb = nan(Nrob,Nrob);
    
    stepVerts = hullVerts(:,:,:,step);
    
    %for every pair of robots
    for i = 1:Nrob
        for j = (i+1):Nrob
            %vertex cloud of hulls
            pairCloud = [stepVerts(:,:,i);stepVerts(:,:,j)];
            
            %svm
            SVM = svmtrain(labels,pairCloud,'-q -t 0');
          
            %hyperplane params
            suppVecs = pairCloud(SVM.sv_indices,:);
            w = SVM.sv_coef' * suppVecs;
            
            normw = norm(w);
            currA = w/normw;
            currb = (-1*SVM.rho)/normw;
            
            stepA(:,i,j) = -1*currA;
            stepb(i,j) = currb;
            stepA(:,j,i) = currA;
            stepb(j,i) = -1*currb;

        end
    end
    
    A(:,:,:,step) = stepA;
    b(:,:,step) = stepb;
end

end
