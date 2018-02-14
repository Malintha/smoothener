clear; close all

%% Robots / Path def

% r1_cyl.above = 1;
% r1_cyl.below = 2;
% r1_cyl.radius = 1;
% r1_cyl.height = r1_cyl.below + r1_cyl.above;

r1_cyl = [1,1,2]; %radius,above,below, seperate r1 from r2

% Paths
r1_path = [0, -1, 0;...
           2,0,0; ...
           2,2,2; ...
           2,5,-1; ...
           1,7,-3; ...
           -2,8,1];

r2_path = [2, -1, -5;...
           4,0,-5; ...
           4,2,-7; ...
           4,5,-6; ...
           3,7,-8; ...
           0,8,-6];
% r2_path = [-5,0,0;...
%            -6,2,1;...
%            -7,3,2;...
%            -6,4,3;...
%            -4,6,2;...
%            -5,8,-1];

nsteps = size(r1_path,1) - 1;



%% compute hyperplanes
%pepare data in expected format
paths = zeros(3,size(r1_path,1),2);
paths(:,:,1) = r1_path';
paths(:,:,2) = r2_path';
types = [1,2];

cylinders = zeros(2,2,3);
dummy_cyl = [0.1,0.1,0.1];
cylinders(1,2,:) = r1_cyl;
cylinders(2,1,:) = dummy_cyl;
cylinders(1,1,:) = dummy_cyl;
cylinders(2,2,:) = dummy_cyl;

%hp computation for each timestep;
[A,b] = robot_hp_waypoints(paths,types,cylinders);

%hp computation (treating waypoints as trajectory samples)
[r1_hull] = swept_cyl_verts(r1_cyl,r1_path);
labels = [ones(size(r1_hull,1),1);-1*ones(size(r2_path,1),1)];
trajcloud = [r1_hull;r2_path];
SVM = svmtrain(labels,trajcloud,'-q -t 0');
%hyperplane params
suppVecs = trajcloud(SVM.sv_indices,:);
trajw = SVM.sv_coef' * suppVecs;
trajA = (trajw/norm(trajw));
trajB = -(SVM.rho)/norm(trajw);


%% Swept Volume hull vertices
[r1_verts] = swept_cyl_verts(r1_cyl,r1_path);


% ~~~~~~~ All Verts ~~~~~~~
% hull vertices
scatter3(r1_verts(:,1),r1_verts(:,2),r1_verts(:,3),'ko')
hold on;
ax = gca;

% ~~~~~~~ Hull at each step ~~~~~~~
for t = 1:nsteps
    DT = delaunayTriangulation(r1_verts(((t-1)*16+1):(t*16),:));
    [K,v] = convexHull(DT);
    trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),'FaceAlpha',0.2,'FaceColor',[0.2,0.5,0.5],'edgecolor','none');
end

% ~~~~~~~ Cylinders at waypoints ~~~~~~~
%   initial cylinder and color scheme
[r1_c1x, r1_c1y, r1_c1z] = robot_cylinder(r1_cyl,r1_path(1,:));
r1_color = zeros(size(r1_c1z),'like',r1_c1z);
r1_color(:,:,1) = 0.1; % red
r1_color(:,:,2) = 0.6; % green
r1_color(:,:,3) = 0.2; % blue

surf(ax,r1_c1x, r1_c1y, r1_c1z,r1_color);
% other cylinders
for p = 2:size(r1_path,1)
    [r1_c2x, r1_c2y, r1_c2z] = robot_cylinder(r1_cyl,r1_path(p,:));
    surf(ax,r1_c2x, r1_c2y, r1_c2z,r1_color);
end

% ~~~~~~~ Paths ~~~~~~~
plot3(r1_path(:,1),r1_path(:,2),r1_path(:,3),'-go', ...
    'LineWidth', 3);
plot3(r2_path(:,1),r2_path(:,2),r2_path(:,3),'-ro', ...
    'LineWidth', 3);


xlabel('x')
ylabel('y')
zlabel('z')
ax.Projection = 'perspective';
ax.DataAspectRatioMode = 'manual';
ax.DataAspectRatio = [1 1 1];
axis vis3d;

% ~~~~~~~ Hyperplanes "Animation" ~~~~~~~
%bounds for hp surf
xrange = [-8,4];
yrange = [-2,9];
zrange = [-3,5];
xyzstep = 1;
%timestep hyperplanes
for t = 1:nsteps
    k = input('next: ');
    if (t>1) %hide plots from last step
        qlast.Visible = 'off';
        slast.Visible = 'off';
        p1last.Visible = 'off';
        p2last.Visible = 'off';
    end
    %plot plane
    Astep = A(:,1,2,t);
    bstep = b(1,2,t);
    [hpx,hpy,hpz] = hyperplane_surf(Astep,bstep,xrange,yrange,zrange,xyzstep);
    u = Astep(1)*ones(size(hpx,1),size(hpx,2));
    v = Astep(2)*ones(size(hpx,1),size(hpx,2));
    w = Astep(3)*ones(size(hpx,1),size(hpx,2));
    qlast = quiver3(ax,hpx,hpy,hpz,u,v,w,0.5);
    slast = surf(ax,hpx,hpy,hpz,'FaceAlpha',0.5,'FaceColor',[0.4,0.1,0.4],'edgecolor','none');
    %highlight segment
    p1last = plot3(r1_path(t:(t+1),1),r1_path(t:(t+1),2),r1_path(t:(t+1),3),'-go', 'LineWidth', 8);
    p2last = plot3(r2_path(t:(t+1),1),r2_path(t:(t+1),2),r2_path(t:(t+1),3),'-ro', 'LineWidth', 8);
end
%full trajectory hyperplane
k = input('Full traj: ');
qlast.Visible = 'off';
slast.Visible = 'off';
p1last.Visible = 'off';
p2last.Visible = 'off';
[hpx,hpy,hpz] = hyperplane_surf(trajA,trajB,xrange,yrange,zrange,xyzstep);
u = trajA(1)*ones(size(hpx,1),size(hpx,2));
v = trajA(2)*ones(size(hpx,1),size(hpx,2));
w = trajA(3)*ones(size(hpx,1),size(hpx,2));
qlast = quiver3(ax,hpx,hpy,hpz,u,v,w,0.5);
slast = surf(ax,hpx,hpy,hpz,'FaceAlpha',0.5,'FaceColor',[0.4,0.1,0.4],'edgecolor','none');

hold off;

%% Helpers
function [x,y,z] = robot_cylinder(rcyl,pos)
    
    %cylinder radius r, height 1 with bottom at origin
    [cx, cy, cz] = cylinder(rcyl(1));
    
    %translate xy, scale and translate z
    x = cx + pos(1);
    y = cy + pos(2);
    z = (cz * (rcyl(2)+rcyl(3))) + pos(3) - rcyl(3);

end

function [x,y,z] = hyperplane_surf(A,b,xrange,yrange,zrange,xyzstep)

    maxdim = find(abs(A)==max(abs(A)));
    %solve for dim that has largest normal coeff (avoid 1/coeff making huge
    %numbers)
    if (maxdim == 1)
        [z, y] = meshgrid(zrange(1):xyzstep:zrange(2),...
                    yrange(1):xyzstep:yrange(2)); % Generate z and y data
        x = (-1/A(1))*(A(3)*z + A(2)*y + b); % Solve for x data
    elseif (maxdim == 2)
        [x, z] = meshgrid(xrange(1):xyzstep:xrange(2),...
                    zrange(1):xyzstep:zrange(2)); % Generate x and z data
        y = (-1/A(2))*(A(1)*x + A(3)*z + b); % Solve for y data
    else
        [x, y] = meshgrid(xrange(1):xyzstep:xrange(2),...
                    yrange(1):xyzstep:yrange(2)); % Generate x and y data
        z = (-1/A(3))*(A(1)*x + A(2)*y + b); % Solve for z data
    end
    
    
    


end