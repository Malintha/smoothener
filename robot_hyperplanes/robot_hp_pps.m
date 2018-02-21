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

%% Compute hyperplanes via SVM

%initialize output structures
A = nan(3,Nrob,Nrob,Nsteps);
b = nan(Nrob,Nrob,Nsteps);

%for each timestep
for step = 1:Nsteps
    %renames for parfor
    stepA = nan(3,Nrob,Nrob);
    stepb = nan(Nrob,Nrob);
    
    %compute svm for every pair of robots
    for i = 1:Nrob
        for j = (i+1):Nrob
            %sample from trajectory
            traj_i = pp_sample_piece(pps{i},step,Neval);
            traj_j = pp_sample_piece(pps{j},step,Neval);
            
            %SHP constraint for j
            %compute conflict hull from i's perspective
            %   j's path must stay out of this hull
            if (step == 4)
                debug = 0;
            end
            [hull] = swept_cyl_verts(cylinders(types(j),types(i),:), traj_i');
                                 
            %vertex cloud for hull + waypoints for agent j
            pairCloud = [hull; traj_j'];
            
            %labels for cloud, 1 for robot i, -1 for j
            labels = [ones(size(hull,1),1);-1*ones(size(traj_j',1),1)];

            %train svm to get hyperplane
            SVM = svmtrain(labels,pairCloud,'-c 10000 -q -t 0');

            %hyperplane params
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
%             if (step ==4)
%                 close all;
%                 scatter3(hull(:,1),hull(:,2),hull(:,3),'bo')
%                 hold on;
%                 plot3(traj_i(1,:),traj_i(2,:),traj_i(3,:),'-go','LineWidth',7)
%                 plot3(traj_j(1,:),traj_j(2,:),traj_j(3,:),'-ro','LineWidth',7)
%                 plotA = -currA;
%                 plotb = currb;
%                 buf = 2;
%                 bbx = [min([traj_i(1,:),traj_j(1,:)])-buf,max([traj_i(1,:),traj_j(1,:)])+buf];
%                 bby = [min([traj_i(2,:),traj_j(2,:)])-buf,max([traj_i(2,:),traj_j(2,:)])+buf];
%                 bbz = [min([traj_i(3,:),traj_j(3,:)])-buf,max([traj_i(3,:),traj_j(3,:)])+buf];
%                 [debx,deby,debz] = hyperplane_surf(plotA,plotb,bbx,bby,bbz,2);
%                 u = plotA(1)*ones(size(debx,1),size(debx,2));
%                 v = plotA(2)*ones(size(deby,1),size(deby,2));
%                 ww = plotA(3)*ones(size(debz,1),size(debz,2));
%                 quiver3(debx,deby,debz,u,v,ww,0.1);
%                 surf(debx,deby,debz,'FaceAlpha',0.5,'FaceColor',[0.4,0.1,0.4],'edgecolor','none');
%                 ax = gca;
%                 xlabel('x')
%                 ylabel('y')
%                 zlabel('z')
%                 ax.Projection = 'perspective';
%                 ax.DataAspectRatioMode = 'manual';
%                 ax.DataAspectRatio = [1 1 1];
%                 axis vis3d;
%                 hold off;
%                 debug = 0;
%             end
            if any(suppLab~=sign(suppDists))
                warning(sprintf('Robots (%d,%d) trajectories conflict at step %d',i,j,step));
                %DEBUG PLOT
                % agent i is green, agent j is red.
                % constraint should be between green/red and pointing
                % towards red
                close all;
                scatter3(hull(:,1),hull(:,2),hull(:,3),'bo')
                hold on;
                plot3(traj_i(1,:),traj_i(2,:),traj_i(1,:),'-go','LineWidth',7)
                plot3(traj_j(1,:),traj_j(2,:),traj_j(1,:),'-ro','LineWidth',7)
                plotA = -currA;
                plotb = currb;
                buf = 2;
                bbx = [min([traj_i(1,:),traj_j(1,:)])-buf,max([traj_i(1,:),traj_j(1,:)])+buf];
                bby = [min([traj_i(2,:),traj_j(2,:)])-buf,max([traj_i(2,:),traj_j(2,:)])+buf];
                bbz = [min([traj_i(3,:),traj_j(3,:)])-buf,max([traj_i(3,:),traj_j(3,:)])+buf];
                [debx,deby,debz] = hyperplane_surf(plotA,plotb,bbx,bby,bbz,2);
                u = plotA(1)*ones(size(debx,1),size(debx,2));
                v = plotA(2)*ones(size(deby,1),size(deby,2));
                ww = plotA(3)*ones(size(debz,1),size(debz,2));
                quiver3(debx,deby,debz,u,v,ww,0.1);
                surf(debx,deby,debz,'FaceAlpha',0.5,'FaceColor',[0.4,0.1,0.4],'edgecolor','none');
                ax = gca;
                xlabel('x')
                ylabel('y')
                zlabel('z')
                ax.Projection = 'perspective';
                ax.DataAspectRatioMode = 'manual';
                ax.DataAspectRatio = [1 1 1];
                axis vis3d;
                hold off;
                debug = 0;
            end
            
            %SHP constraint for i
            %compute conflict hull from j's perspective
            %   i's path must stay out of this hull
            [hull] = swept_cyl_verts(cylinders(types(i),types(j),:),traj_j');
                                 
            %vertex cloud for hull + waypoints for agent i
            pairCloud = [hull; traj_i'];
            
            %labels for cloud, 1 for robot j, -1 for i
            labels = [ones(size(hull,1),1);-1*ones(size(traj_i',1),1)];

            %train svm to get hyperplane
            SVM = svmtrain(labels,pairCloud,'-c 1000 -q -t 0');

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
                warning(sprintf('Robots (%d,%d) trajectories conflict at step %d',i,j,step));
            end

        end
    end
    
    A(:,:,:,step) = stepA;
    b(:,:,step) = stepb;
end


end

