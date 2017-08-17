function [ A, b ] = robot_hp_pps( pps, ellipsoids, Neval )
%ROBOT_HP_PPS Calculates robot seperating hyperplanes
%INPUT:
%   pps: {#robots} cell array of Matlab ppform structs
%         pps{n} is the ppform for agent n
%   ellipsoids: [#robots, #dimensions]
%          ellipsoids(n,:) the axis-aligned ellipsoid axes size for
%          robot n
%   neval : integer specifying number of samples for each piece of pps for
%           generating swept hull
%OUTPUT:
%   A: [DIM x NROBOTS x NROBOTS x (NPTS - 1)] array of 
%      hyperplane normal vectors for each robot-robot interaction
%      at each segment.
%   b: distance from origin for hyperplanes. i.e. a(:,...)^T x <= b
%NOTES:
%   currently only works for 3d paths / ellipsoids


%% 

Nrob = length(pps); % #robots
Nsteps = pps{1}.pieces; % #time steps
dim = pps{1}.dim;
Nbreaks = pps{1}.breaks;

%% generate untransformed-ellipsoid vertices
ntheta = 11;
nphi = 11;
%unit sphere
sphereVerts = ellipsoid_verts(1,1,1,ntheta,nphi);
%get ellipsoids by scaling sphere
ellVerts = zeros(size(sphereVerts,1),3,Nrob);
for r = 1:Nrob
    ellVerts(:,:,r) = [sphereVerts(:,1).*ellipsoids(r,1),...
                       sphereVerts(:,2).*ellipsoids(r,2),...
                       sphereVerts(:,3).*ellipsoids(r,3)];
end

 %TODO ABOVE: uniform meshes for ellipsoids. 
 %  http://persson.berkeley.edu/distmesh/

NellVerts = size(ellVerts,1); %number of vertices per ellipsoid

%% generate trajectory hulls
%   by transforming ellipsoids along sampled pps

for i=2:Nrob
    assert(pps{i}.dim == dim);
    assert(pps{i}.pieces == Nsteps);
    assert(all(abs(pps{i}.breaks - Nbreaks) < 10e-9));
end

%translate ellipsoids to create 'swept hull' vertices for each timestep
hullVerts = zeros(Neval*NellVerts,size(ellVerts,2),Nrob, Nsteps);
for r = 1:Nrob
    for step = 1:Nsteps
        %sample pps
        transforms = pp_sample_piece(pps{r},step,Neval);
        %for each sample, place ellipsoid
        for t = 1:Neval
            tf = repmat(transforms(:,t)',NellVerts,1);
            s = (t-1)*NellVerts + 1; 
            e = t*NellVerts;
            hullVerts(s:e,:,r,step) = ellVerts(:,:,r) + tf;
        end     
    end
end

%% Compute hyperplanes via SVM

%labels for svm
%   true for robot 1, false for robot 2
labels = [true(Neval*NellVerts,1);false(Neval*NellVerts,1)];

%initialize output structures
A = nan(3,Nrob,Nrob,Nsteps);
b = nan(Nrob,Nrob,Nsteps);

%for each timestep
parfor step = 1:Nsteps
    %renames for parfor
    stepA = nan(3,Nrob,Nrob);
    stepb = nan(Nrob,Nrob);
    stepVerts = hullVerts(:,:,:,step);
    
    %compute svm for every pair of robots
    for i = 1:Nrob
        for j = (i+1):Nrob
            %svm linear classifier
            SVM = fitcsvm([stepVerts(:,:,i);stepVerts(:,:,j)], labels);

            %hyperplane params
            currA = SVM.Beta ./ norm(SVM.Beta);
            currb = SVM.Bias ./ norm(SVM.Beta);
            
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

