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
            [hull] = swept_cyl_verts(cylinders(types(i),types(j),:),...
                                     [paths(:,step,i)';paths(:,step+1,i)']);
                                 
            %vertex cloud for hull + waypoints for agent j
            pairCloud = [hull; paths(:,step,j)';paths(:,step+1,j)'];
            
            %labels for cloud, 1 for robot i, -1 for j
            labels = [ones(size(hull,1),1);-1;-1];
            
            %DEEEEEEEEEEEEEEBBBBBBUUUUUUUUUUUGGGGGGGGGGGg
%             pairCloud = [paths(:,step,i)';paths(:,step+1,i)'; paths(:,step,j)';paths(:,step+1,j)'];
%             labels = [1;1;-1;-1];
            
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
            
            %DEEEEEEEEEEEEEBBBBBBBBBBBBUUUUUUUUUUUGGGGGGGGGGGGG
%             close all
%             scatter3(hull(:,1),hull(:,2),hull(:,3),'go')
%             hold on;
%             plot3([paths(1,step,i);paths(1,step+1,i)],[paths(2,step,i)';paths(2,step+1,i)],[paths(3,step,i)';paths(3,step+1,i)],'-go')
%             plot3([paths(1,step,j);paths(1,step+1,j)],[paths(2,step,j)';paths(2,step+1,j)],[paths(3,step,j)';paths(3,step+1,j)],'-ro')
%             plotA = -currA;
%             plotb = currb;
%             [debx,deby,debz] = hyperplane_surf(plotA,plotb,[-5,5],[-5,5],[-1,3],2);
%             u = plotA(1)*ones(size(debx,1),size(debx,2));
%             v = plotA(2)*ones(size(deby,1),size(deby,2));
%             ww = plotA(3)*ones(size(debz,1),size(debz,2));
%             quiver3(debx,deby,debz,u,v,ww,0.1);
%             surf(debx,deby,debz,'FaceAlpha',0.5,'FaceColor',[0.5,0.1,0.4],'edgecolor','none');
%             debug = 0;
            
            %SHP constraint for i
            %compute conflict hull from j's perspective
            %   i's path must stay out of this hull
            [hull] = swept_cyl_verts(cylinders(types(j),types(i),:),...
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
            
            %DEEEEEEEEEEEEEBBBBBBBBBBBBUUUUUUUUUUUGGGGGGGGGGGGG
%             scatter3(hull(:,1),hull(:,2),hull(:,3),'ro')
%             plotA = -currA;
%             plotb = currb;
%             [debx,deby,debz] = hyperplane_surf(plotA,plotb,[-5,5],[-5,5],[-1,3],2);
%             u = plotA(1)*ones(size(debx,1),size(debx,2));
%             v = plotA(2)*ones(size(deby,1),size(deby,2));
%             ww = plotA(3)*ones(size(debz,1),size(debz,2));
%             quiver3(debx,deby,debz,u,v,ww,0.1);
%             surf(debx,deby,debz,'FaceAlpha',0.5,'FaceColor',[0.1,0.5,0.4],'edgecolor','none');
%             ax = gca;
%             xlabel('x')
%             ylabel('y')
%             zlabel('z')
%             ax.Projection = 'perspective';
%             ax.DataAspectRatioMode = 'manual';
%             ax.DataAspectRatio = [1 1 1];
%             axis vis3d;
%             debug = 0;
            
        end
    end
    
    A(:,:,:,step) = stepA;
    b(:,:,step) = stepb;
end

end