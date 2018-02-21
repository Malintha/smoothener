function [ A, b ] = robot_hp_waypoints( paths, types,cylinders )
%ROBOT_HP_WAYPOINTS Calculates robot seperating hyperplanes
%INPUT:
%   paths: [3, #waypoints, #robots]
%       paths(:,:,n) are the waypoints for robot n
%   types: [n]
%       types(n) is the type of robot traversing paths(:,:,n).
%       Integer in range [1,#types]
%   cylinders: [#types,#types, 3]
%       cylinders(a,b,1) radius of cyl needed to separate type 'a' from 'b'
%       cylinders(a,b,2) height 'a' needs to be above 'b'
%       cylinders(a,b,3) height 'a' needs to be below 'b'
%OUTPUT:
%   A: [DIM x NROBOTS x NROBOTS x (NPTS - 1)] array of 
%      hyperplane normal vectors for each robot-robot interaction
%      at each segment.
%   b: distance from origin for hyperplanes. i.e. a(:,...)^T x <= b
%NOTES:
%   currently only works for 3d paths / ellipsoids

Nrob = size(paths,3); % #robots
Nsteps = size(paths,2) -1; % #time steps

%initialize output structures
A = nan(3,Nrob,Nrob,Nsteps);
b = nan(Nrob,Nrob,Nsteps);

%parfor each timestep
for step = 1:Nsteps
    %private renames for parfor
    stepA = nan(3,Nrob,Nrob);
    stepb = nan(Nrob,Nrob);
    
    %stepVerts = hullVerts(:,:,:,step);
    %for every pair of robots
    for i = 1:Nrob
        for j = (i+1):Nrob
            %SHP constraint for j
            %compute conflict hull from i's perspective
            %   j's path must stay out of this hull
            [hull] = swept_cyl_verts(cylinders(types(j),types(i),:),...
                                     [paths(:,step,i)';paths(:,step+1,i)']);
                                 
            %vertex cloud for hull + waypoints for agent j
            pairCloud = [hull; paths(:,step,j)';paths(:,step+1,j)'];
            
            %labels for cloud, 1 for robot i, -1 for j
            labels = [ones(size(hull,1),1);-1;-1];
            
            %train svm to get hyperplane
            SVM = svmtrain(labels,pairCloud,'-c 10000 -q -t 0');
            %extract params
            suppVecs = pairCloud(SVM.sv_indices,:);
            w = SVM.sv_coef' * suppVecs;
            normw = norm(w);
            currA = w/normw;
            currb = (SVM.rho)/normw;
            
            stepA(:,j,i) = currA;
            stepb(j,i) = currb;
            
            %sanity check for inseperable case
            suppDists = suppVecs * currA' - currb;
            suppLab = labels(SVM.sv_indices,:);
            if any(suppLab~=sign(suppDists))
                warning(sprintf('Robots (%d,%d) paths conflict at step %d',i,j,step));
            end
            
            %DEEEEEEEEEEEEEBBBBBBBBBBBBUUUUUUUUUUUGGGGGGGGGGGGG
%             traj_i = [paths(:,step,i)';paths(:,step+1,i)']';
%             traj_j = [paths(:,step,j)';paths(:,step+1,j)']';
%             % agent i is green, agent j is red.
%             % constraint should be between green/red and pointing
%             % towards red
%             close all;
%             scatter3(hull(:,1),hull(:,2),hull(:,3),'bo')
%             hold on;
%             plot3(traj_i(1,:),traj_i(2,:),traj_i(3,:),'-go','LineWidth',7)
%             plot3(traj_j(1,:),traj_j(2,:),traj_j(3,:),'-ro','LineWidth',7)
%             plotA = -currA;
%             plotb = currb;
%             buf = 2;
%             bbx = [min([traj_i(1,:),traj_j(1,:)])-buf,max([traj_i(1,:),traj_j(1,:)])+buf];
%             bby = [min([traj_i(2,:),traj_j(2,:)])-buf,max([traj_i(2,:),traj_j(2,:)])+buf];
%             bbz = [min([traj_i(3,:),traj_j(3,:)])-buf,max([traj_i(3,:),traj_j(3,:)])+buf];
%             [debx,deby,debz] = hyperplane_surf(plotA,plotb,bbx,bby,bbz,2);
%             u = plotA(1)*ones(size(debx,1),size(debx,2));
%             v = plotA(2)*ones(size(deby,1),size(deby,2));
%             ww = plotA(3)*ones(size(debz,1),size(debz,2));
%             quiver3(debx,deby,debz,u,v,ww,0.1);
%             surf(debx,deby,debz,'FaceAlpha',0.5,'FaceColor',[0.4,0.1,0.4],'edgecolor','none');
%             ax = gca;
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             ax.Projection = 'perspective';
%             ax.DataAspectRatioMode = 'manual';
%             ax.DataAspectRatio = [1 1 1];
%             axis vis3d;
%             hold off;
%             debug = 0;
            
            %SHP constraint for i
            %compute conflict hull from j's perspective
            %   i's path must stay out of this hull
            [hull] = swept_cyl_verts(cylinders(types(i),types(j),:),...
                                     [paths(:,step,j)';paths(:,step+1,j)']);

            %vertex cloud for hull + waypoints for agent i
            pairCloud = [hull; paths(:,step,i)';paths(:,step+1,i)'];
            
            %labels for cloud, 1 for robot j, -1 for i
            labels = [ones(size(hull,1),1);-1;-1];

            %train svm to get hyperplane
            SVM = svmtrain(labels,pairCloud,'-c 10000 -q -t 0');

            %hyperplane params
            suppVecs = pairCloud(SVM.sv_indices,:);
            w = SVM.sv_coef' * suppVecs;
            normw = norm(w);
            currA = w/normw;
            currb = (SVM.rho)/normw;
            stepA(:,i,j) = currA;
            stepb(i,j) = currb;
            
            %sanity check for inseperable case
            suppDists = suppVecs * currA' - currb;
            suppLab = labels(SVM.sv_indices,:);
            if any(suppLab~=sign(suppDists))
                warning(sprintf('Robots (%d,%d) paths conflict at step %d',i,j,step));
            end
            
            %DEEEEEEEEEEEEEBBBBBBBBBBBBUUUUUUUUUUUGGGGGGGGGGGGG
            
        end
    end
    
    A(:,:,:,step) = stepA;
    b(:,:,step) = stepb;
end

end