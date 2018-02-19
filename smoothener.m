%% Problem Specification

% Need 10 inputs. Definitions
%   k = number of waypoints in discrete plan
%   N = number of agents
%   T = number of agent types

%INPUTS: 
% paths: [dim k N] float
%   paths(:,:,i) k waypoints for agent i

% types: [N] int
%   types(i) number in {1...T} for agent i 

% conf_cylinders: [T,T,3] float
%   defines conflict geometry for every pair in T types
%   (i,j,1) separating radius of type i and type j
%   (i,j,2) type i must be this far above type j
%   (i,j,3) type i must be this far below type j

% obs_cylinders: [T,3] float
%   defines conflict geometry for every pair of T types
%   (i,j,1) separating radius of type i and obstacles
%   (i,j,2) agent i must be this far above obstacle
%   (i,j,3) agent i must be this far below obstacle

% deg: [1] int
%   polynomial degree

% cont: [1] int
%   derivative continuity

% timescale: [1] float
%   seconds per timestep in discrete plan

% iters: [1] int
%   number of smoothing iterations

% pp_obs_sep_fun: [1] func handle
%   obstacle separating function

%% Input Specification
clear; close all;
% r1_path = [0, 0,  3;...
%            0, 2,  3; ...
%            0, 4,  3];
% 
% r2_path = [4, 0,  3;...
%            4, 2,  3;...
%            4, 4,  3];
% 
% r3_path = [-4, 0,  3;...
%            -4, 2,  3; ...
%            -4, 4,  3];
% 
% r4_path = [0, 0,  7;...
%            0, 2,  7;...
%            0, 4,  7];
%        
% r5_path = [0, 0,  0;...
%            0, 2,  0; ...
%            0, 4,  0];

% ~~~~~~ Env file for octomap ~~~~~~
map = '/home/mark/act/smoothener/examples/swap4/map.bt';

% ~~~~~~ paths Input ~~~~~~
paths = read_schedule('./examples/swap4/discreteSchedule.json');
[dim, k, N] = size(paths);
nsteps = size(paths,2)-1;

% ~~~~~~ types Input ~~~~~~
%1 = small, 2 = large for swap4
ntypes = 2;
types = [1;1;2;2];

% types = [1;2];
% ~~~~~~ conf_cylinders Input ~~~~~~
conf_cylinders = zeros(ntypes,ntypes,3);
%cylinders(i,j,1) = radius type i must stay away from type j
%cylinders(i,j,2) = radius type i must stay above type j
%cylinders(i,j,3) = radius type i must stay below type j

conf_cylinders(1,2,:) = [0.20,0.30,0.60];
conf_cylinders(2,1,:) = [0.20,0.60,0.30];

conf_cylinders(1,1,:) = [0.15,0.30,0.30];
conf_cylinders(2,2,:) = [0.25,0.50,0.50];

% ~~~~~~ obs_cylinders Input NOTE: CURRENTLY USING ELLIPSOIDS ~~~~~~
obs_cylinders = ones(ntypes,3);
%obs_cylinders(i,:) = [radius,above,below] for environment
%Right now it is [rx,ry,rz] for ellipsoids
obs_cylinders(1,:) = [0.15,0.15,0.15]*.1;
obs_cylinders(2,:) = [0.25,0.25,0.25]*.1;
%hack for ellipsoid input to octomap separation function
obs_ellipsoids = zeros(N,3);
for n = 1:N
    obs_ellipsoids(n,:) = obs_cylinders(types(n),:);
end

% ~~~~~~ bbox Input ~~~~~~
bbox = [-5.5, 2.0;...
        -3.0, 3.5;...
         0.0, 2.5];
     
% bbbuffer = 10;
% bbox = [min(min(paths(1,:,:))) - bbbuffer,max(max(paths(1,:,:)) + bbbuffer);...
%         min(min(paths(2,:,:))) - bbbuffer,max(max(paths(2,:,:)) + bbbuffer);...
%         min(min(paths(3,:,:))) - bbbuffer,max(max(paths(3,:,:)) + bbbuffer)];
    
% ~~~~~~ deg,cont,timescale,iters input ~~~~~~
deg = 7;
cont = 4;
timescale = 1;
iters = 2;
Neval = 32; %number of samples on pps separation

% ~~~~~~ pp obstacle separation function ~~~~~~
pp_obs_sep_fun = @(poly,elip) pp_obs_sep_octomap(poly,elip,map);

%% Smoothener
[dim, k, N] = size(paths);

% for a reasonable problem, cost should converge after ~5 iterations.
assert(iters >= 1);
assert(iters <= 20);

% outputs
all_pps = cell(iters, N);
all_costs = zeros(iters, N);
all_corridors = cell(iters, N);

% piecewise linear (physically impossible) pps of path
% for input to pp-vs-octree obstacle hyperplane function
pps = path_linear_pps(paths, timescale, deg + 1);

for iter=1:iters
    fprintf('iteration %d of %d...\n', iter, iters);
    tic;
    if iter==1
        % first iteration: decompose by segments
        [A, b] = robot_hp_waypoints(paths, types, conf_cylinders);
    else
        % continuing iteration: decompose by pps
        [A, b] = robot_hp_pps(pps,types,conf_cylinders,Neval);
    end

    if iter > 1
        for irobot=1:N
            paths(:,:,irobot) = ppval(pps{irobot}, pps{irobot}.breaks);
        end
    end

    hs = pp_obs_sep_fun(pps, obs_ellipsoids);

    t_hyperplanes = toc;
    fprintf('hyperplanes: %f sec\n', t_hyperplanes);

    % solve the independent spline trajectory optimization problems.
    tic;
    pps = cell(1,N);
    iter_costs = zeros(1,N);
    % parfor
    for j=1:N
        fprintf(' agent %d of %d...\n', j, N);
        lb = bbox(:,1) + [obs_cylinders(types(j),1);obs_cylinders(types(j),1);obs_cylinders(types(j),3)];
        ub = bbox(:,2) - [obs_cylinders(types(j),1);obs_cylinders(types(j),1);obs_cylinders(types(j),2)];

        hs_slice = squeeze(hs(j,:));
        step_n_faces = cellfun(@(a) size(a, 1), hs_slice);
        assert(length(step_n_faces) == (k-1));
        max_n_faces = max(step_n_faces(:));

        Aobs = nan(dim, max_n_faces, k-1);
        bobs = nan(max_n_faces, k-1);
        Arobots = squeeze(A(:,j,:,:));
        brobots = squeeze(b(j,:,:));
        for i=1:(k-1)
            n_faces = step_n_faces(i);
            hs_slice_step = hs_slice{i};
            assert(size(hs_slice_step, 2) == 4);
            Aobs(:,1:n_faces,i) = hs_slice_step(:,1:3)';
            bobs(1:n_faces,i) = hs_slice_step(:,4);
        end
        [pps{j}, iter_costs(j)] = corridor_trajectory_optimize(...
            Arobots, brobots, ...
            Aobs, bobs, ...
            lb, ub,...
            paths(:,:,j), deg, cont, timescale, obs_cylinders(types(j),:));%[0.2 0.2 0.4], [0.2 0.2 0.2]);%

        s = [];
        s.Arobots = Arobots;
        s.brobots = brobots;
        s.Aobs = Aobs;
        s.bobs = bobs;
        all_corridors{iter,j} = s;
    end
    t_splines = toc;
    fprintf('splines: %f sec\n', t_splines);
    fprintf('total: %f sec\n', t_hyperplanes + t_splines);
    fprintf('cost: %f\n', sum(iter_costs));
    all_costs(iter,:) = iter_costs;
    all_pps(iter,:) = pps;
end

%% Plot?

for i=1:size(paths,3)	% [14,16,17,31]		
    h = plot3n(paths(:,:,i));
    hold on;
    color = get(h, 'color');
    %delete(h);
    duration = pps{i}.breaks(end);
    t = 0:0.05:duration;
    x = ppval(pps{i}, t);
    h = plot3n(x, 'color', color, 'LineWidth', 3);
%     lastPos = ppval(pps{i}, duration);
%     scatter3(lastPos(1), lastPos(2), lastPos(3),50,'k','filled');
%     firstPos = ppval(pps{i}, 0);
%     scatter3(firstPos(1), firstPos(2), firstPos(3),50,'k', 'square');%'g','filled');
end
