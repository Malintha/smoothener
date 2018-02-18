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

%Path defs for 2 robots
% r1_path = [4, 1,  0;...
%            3, 2,  1; ...
%            2, 3,  2; ...
%            1, 4,  3; ...
%            0, 5,  4; ...
%            -1, 6,  5];
% 
% r2_path = [-4, 1,  0;...
%            -3, 2,  -1;...
%            -2, 3,  -2;...
%            -1, 4,  -3;...
%            -0, 5,  -4;...
%            1, 6,  -5];


r1_path = [0, 0,  3;...
           0, 2,  3; ...
           0, 4,  3];

r2_path = [4, 0,  3;...
           4, 2,  3;...
           4, 4,  3];
% r1_path = [0, 0,  3;...
%            0, 2,  3]
% 
% r2_path = [4, 0,  3;...
%            4, 2,  3]
r3_path = [-4, 0,  3;...
           -4, 2,  3; ...
           -4, 4,  3];

r4_path = [0, 0,  7;...
           0, 2,  7;...
           0, 4,  7];
       
r5_path = [0, 0,  0;...
           0, 2,  0; ...
           0, 4,  0];

       
%pad
% r1_path = [r1_path(1,:);r1_path;r1_path(3,:)]
% r2_path = [r2_path(1,:);r2_path;r2_path(3,:)]
       
assert(size(r1_path,1) == size(r2_path,1))

nsteps = size(r1_path,1) - 1;

% ~~~~~~ paths Input ~~~~~~
paths = zeros(3,size(r1_path,1),5);
paths(:,:,1) = r1_path';
paths(:,:,2) = r2_path';
paths(:,:,3) = r3_path';
paths(:,:,4) = r4_path';
paths(:,:,5) = r5_path';

% ~~~~~~ types Input ~~~~~~
types = [1;2;2;2;2];
% types = [1;2];
% ~~~~~~ conf_cylinders Input ~~~~~~
conf_cylinders = zeros(2,2,3);
% currently only separating type 1 from 2
conf_cylinders(1,2,:) = [1,1,1];
conf_cylinders(2,1,:) = [1,1,1];
conf_cylinders(1,1,:) = [1,1,1];
conf_cylinders(2,2,:) = [1,1,1];

% ~~~~~~ obs_cylinders Input NOTE: CURRENTLY USING ELLIPSOIDS ~~~~~~
obs_cylinders = ones(5,3);
obs_ellipsoids = obs_cylinders;

% ~~~~~~ bbox Input ~~~~~~
bbbuffer = 10;
bbox = [min(min(paths(1,:,:))) - bbbuffer,max(max(paths(1,:,:)) + bbbuffer);...
        min(min(paths(2,:,:))) - bbbuffer,max(max(paths(2,:,:)) + bbbuffer);...
        min(min(paths(3,:,:))) - bbbuffer,max(max(paths(3,:,:)) + bbbuffer)];
    
% ~~~~~~ deg,cont,timescale,iters input ~~~~~~
deg = 7;
cont = 4;
timescale = 1;
iters = 2;
Neval = 32; %number of samples on pps separation

% ~~~~~~ pp obstacle separation function ~~~~~~
pp_obs_sep_fun = @pp_obs_sep_none;

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
        lb = bbox(:,1) + [obs_cylinders(j,1);obs_cylinders(j,1);obs_cylinders(j,3)];
        ub = bbox(:,2) - [obs_cylinders(j,1);obs_cylinders(j,1);obs_cylinders(j,2)];

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
        [pps{j}, iter_costs(j)] = corridor_trajectory_optimize_devel(...
            Arobots, brobots, ...
            Aobs, bobs, ...
            lb, ub,...
            paths(:,:,j), deg, cont, timescale, obs_cylinders(j,:));%[0.2 0.2 0.4], [0.2 0.2 0.2]);%

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
