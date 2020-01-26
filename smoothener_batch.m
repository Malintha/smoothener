% Function for batch processing
%
% inputs:
%   map: env file for octomap (bt)
%   schedule_file: file for discrete schedule (yaml)
%   types file: file with robot types (yaml)
%   outcsv: output folder
% outputs:
%   None (writes csv files in output folder)
%
function smoothener_batch(map, schedule_file, types_file, outcsv)

    % ~~~~~~ USER INPUT ~~~~~~
    OBSTACLES = true;

    % ~~~~~~ deg,cont,timescale,iters input ~~~~~~
    deg = 7;
    cont = 1;%4;
    timescale = 1.0;
    iters = 8;
    Neval = 32; %number of samples on pps separation

    % ~~~~~~ END USER INPUT ~~~~~~~~~~~~

    [paths,names,typeNames] = read_schedule(schedule_file);
    [dim, k, N] = size(paths);
    nsteps = size(paths,2)-1;

    % plot3(paths(1,:,1),paths(2,:,1),paths(3,:,1))
    % hold on
    % plot3(paths(1,:,2),paths(2,:,2),paths(3,:,2))

    % ~~~~~ read types ~~~~~~
    typesStruct = yaml.ReadYaml(types_file);

    ntypes = size(typesStruct.agentTypes);
    ntypes = ntypes(2);
    locomotion = ones(ntypes, 1) * 3;

    %Right now it is [rx,ry,rz] for ellipsoids
    obs_cylinders = ones(ntypes,3);

    % fill obs_cylinders and locomotion
    for i=1:ntypes
        shape = typesStruct.agentTypes{1,i}.shape;
        if strcmp(shape.type, "cylinder")
            obs_cylinders(i,:) = [shape.radius, shape.radius, shape.height / 2];
        end
        if strcmp(shape.type, "sphere")
            obs_cylinders(i,:) = [shape.radius, shape.radius, shape.radius];
        end
        if contains(typesStruct.agentTypes{1,i}.type, "ground")
            locomotion(i) = 2;
        end
    end

    %obs_cylinders = obs_cylinders * 0.5;

    % fill conf_cylinders
    conf_cylinders = zeros(ntypes,ntypes,3);
    for i=1:ntypes
        for j=1:ntypes
            type_i = typesStruct.agentTypes{1,i}.type;
            type_j = typesStruct.agentTypes{1,j}.type;
            for k=1:length(typesStruct.agentInteractions)
                interaction=typesStruct.agentInteractions{k};
                if strcmp(interaction.typeA, type_i) && strcmp(interaction.typeB, type_j)
                    conf_cylinders(i,j,:) = [interaction.radius, interaction.below, interaction.above];
                end
                if strcmp(interaction.typeA, type_j) && strcmp(interaction.typeB, type_i)
                    conf_cylinders(i,j,:) = [interaction.radius, interaction.above, interaction.below];
                end
            end
        end
    end

    % fill types
    types = zeros(N, 1);

    for n = 1:N
       for i=1:ntypes
         type = typesStruct.agentTypes{1,i}.type;
         if strcmp(typeNames{n}, type)
            types(n) = i;
         end
       end
    end

    %hack for ellipsoid input to octomap separation function
    obs_ellipsoids = zeros(N,3);
    for n = 1:N
        obs_ellipsoids(n,:) = obs_cylinders(types(n),:);
    end

    % ~~~~~~ bbox Input ~~~~~~
    bbox = read_octomap_bbox_mex(map);
    % bbox = [-5.5, 2.0;...
    %         -3.0, 3.5;...
    %          0.0, 2.5];
         
    % bbbuffer = 10;
    % bbox = [min(min(paths(1,:,:))) - bbbuffer,max(max(paths(1,:,:)) + bbbuffer);...
    %         min(min(paths(2,:,:))) - bbbuffer,max(max(paths(2,:,:)) + bbbuffer);...
    %         min(min(paths(3,:,:))) - bbbuffer,max(max(paths(3,:,:)) + bbbuffer)];
        
    % ~~~~~~ pp obstacle separation function ~~~~~~
    if (OBSTACLES)
        pp_obs_sep_fun = @(poly,elip) pp_obs_sep_octomap(poly,elip,map);
    else
        pp_obs_sep_fun = @pp_obs_sep_none;
    end


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
        parfor j=1:N
            fprintf(' agent %d of %d...\n', j, N);
            lb = bbox(:,1) + [obs_cylinders(types(j),1);obs_cylinders(types(j),1);obs_cylinders(types(j),3)];
            ub = bbox(:,2) - [obs_cylinders(types(j),1);obs_cylinders(types(j),1);obs_cylinders(types(j),2)];

            hs_slice = squeeze(hs(j,:));
            step_n_faces = cellfun(@(a) size(a, 1), hs_slice);
            assert(length(step_n_faces) == (k-1));
            max_n_faces = max(step_n_faces(:));

            Aobs = nan(dim, max_n_faces, k-1);
            bobs = nan(max_n_faces, k-1);
            if N == 1
                Arobots = A(:,j,:,:); %squeeze(A(:,j,:,:));
                brobots = b(j,:,:);%squeeze(b(j,:,:));
            else
                Arobots = squeeze(A(:,j,:,:));
                brobots = squeeze(b(j,:,:));
            end
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
                paths(:,:,j), deg, cont, timescale, obs_cylinders(types(j),:),j,iter,locomotion(types(j)));%[0.2 0.2 0.4], [0.2 0.2 0.2]);%
            s = [];
            s.Arobots = Arobots;
            s.brobots = brobots;
            s.Aobs = Aobs;
            s.bobs = bobs;
            all_corridors{iter,j} = s;
        end
        
        if isnan(sum(iter_costs))
            if iter == 1
                display('Failure on first iteration! Exiting...');
                exit;
            end
            pps = all_pps(iter-1,:);
            break;
        end
        
        t_splines = toc;
        fprintf('splines: %f sec\n', t_splines);
        fprintf('total: %f sec\n', t_hyperplanes + t_splines);
        fprintf('cost: %f\n', sum(iter_costs));
        all_costs(iter,:) = iter_costs;
        all_pps(iter,:) = pps;
    end

    %% Save

    % ~~~~~~ pps Output ~~~~~~
    for n = 1:N
        pp2csv(pps{n}, [outcsv,names{n},'.csv'])
    end
end
