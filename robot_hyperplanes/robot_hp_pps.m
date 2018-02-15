function [ A, b ] = robot_hp_pps( pps, types,cylinders, Neval )
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

for i=2:Nrob
    assert(pps{i}.dim == dim);
    assert(pps{i}.pieces == Nsteps);
    assert(all(abs(pps{i}.breaks - Nbreaks) < 10e-9));
end

%translate ellipsoids to create 'swept hull' vertices for each timestep
% hullVerts = zeros(Neval*NellVerts,size(ellVerts,2),Nrob, Nsteps);
% for r = 1:Nrob
%     for step = 1:Nsteps
%         %sample pps
%         transforms = pp_sample_piece(pps{r},step,Neval);
%         %for each sample, place ellipsoid
%         for t = 1:Neval
%             tf = repmat(transforms(:,t)',NellVerts,1);
%             s = (t-1)*NellVerts + 1; 
%             e = t*NellVerts;
%             hullVerts(s:e,:,r,step) = ellVerts(:,:,r) + tf;
%         end     
%     end
% end

%% Compute hyperplanes via SVM

%labels for svm
%   1 for robot 1, -1 for robot 2
labels = [ones(Neval*NellVerts,1);-1*ones(Neval*NellVerts,1)];

%initialize output structures
A = nan(3,Nrob,Nrob,Nsteps);
b = nan(Nrob,Nrob,Nsteps);

%for each timestep
parfor step = 1:Nsteps
    %renames for parfor
    stepA = nan(3,Nrob,Nrob);
    stepb = nan(Nrob,Nrob);
    
    %compute svm for every pair of robots
    for i = 1:Nrob
        for j = (i+1):Nrob
            %sample from trajectory
            traj_i = pp_sample_piece(pps{i},step,Neval);
            traj_j = pp_sample_piece(pps{j},step,Neval);
            
            %compute conflict hull from perspective of agent i
            [hull] = swept_cyl_verts(cylinders(types(i),types(j),:), traj_i');

            %vertex cloud for hull + waypoints for agent j
            pairCloud = [hull; traj_j'];

            %labels for cloud, 1 for robot i, -1 for j
            labels = [ones(size(hull,1),1);-1*ones(size(traj_j,2),1)];

            %train svm to get hyperplane
            SVM = svmtrain(labels,pairCloud,'-q -t 0');

            %hyperplane params
            suppVecs = pairCloud(SVM.sv_indices,:);
            w = SVM.sv_coef' * suppVecs;
            normw = norm(w);
            currA = w/normw;
            currb = (-1*SVM.rho)/normw;

            stepA(:,i,j) = currA;
            stepb(i,j) = currb;
            stepA(:,j,i) = -1*currA;
            stepb(j,i) = -1*currb;

        end
    end
    
    A(:,:,:,step) = stepA;
    b(:,:,step) = stepb;
end


end

