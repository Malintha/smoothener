%% file name
probname = 'act_25';

%% type definitions
types = cell(3,3); %[#types, 3]

types{1,1} = 'cf'; %type name
types{1,2} = [0.15 0.15 0.3]; %conflict size
types{1,3} = [0.15 0.15 0.15]; %obstacle size

types{2,1} = 'mf'; %type name
types{2,2} = [0.2 0.2 0.4]; %conflict size
types{2,3} = [0.2 0.2 0.2]; %obstacle size

types{3,1} = 'bf'; %type name
types{3,2} = [0.25 0.25 0.5]; %conflict size
types{3,3} = [0.25 0.25 0.25]; %obstacle size

%% agent types specification

tid = [3,3,3,3,3,3,3,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2];

%for saving later
agtypes = types(tid,1);


Nrob = numel(agtypes);

%% start / end verts

load('add.mat');
starts = startPos;
goals = goalPos;

%% write types file
fileId = fopen([probname '_types.yaml'],'w');
fprintf(fileId, 'agentTypes:\n');
for i = 1:size(types,1)
    fprintf(fileId,'  - type: %s\n',types{i,1});
    fprintf(fileId,'    conflictSize: [%.2f,%.2f,%.2f]\n',types{i,2}(1),types{i,2}(2),types{i,2}(3));
    fprintf(fileId,'    obstacleSize: [%.2f,%.2f,%.2f]\n',types{i,3}(1),types{i,3}(2),types{i,3}(3));
end
fclose(fileId);

%% write agents file
fileId = fopen([probname '_agents.yaml'],'w');
fprintf(fileId, 'agents:\n');
for i = 1:Nrob
    n = [agtypes{i}, '_', num2str((i > Nrob/2) + 1)];
    fprintf(fileId,'  - name: %s\n',n);
    fprintf(fileId,'    type: %s\n', agtypes{i});
    fprintf(fileId,'    start: start%i\n',i);
    fprintf(fileId,'    goal: goal%i\n',i);
end
fclose(fileId);

%% write add verts

fileId = fopen([probname '_add.yaml'],'w');
fprintf(fileId, 'vertices:\n');

for i = 1:Nrob
    fprintf(fileId,'   - name: start%i\n',i);
    fprintf(fileId,'     pos: [%.3f,%.3f,%.3f]\n',starts(i,1),starts(i,2),starts(i,3));
    fprintf(fileId,'   - name: goal%i\n',i);
    fprintf(fileId,'     pos: [%.3f,%.3f,%.3f]\n',goals(i,1),goals(i,2),goals(i,3));
end
fclose(fileId);

%% figure
close all;
figure();
hold on;
safety = 2;
for i = 1:Nrob
    ellipsoid(starts(i,1),starts(i,2),starts(i,3), ...
      types{tid(i),2}(1)*safety, types{tid(i),2}(2)*safety,types{tid(i),2}(3)*safety);
end
scatter3(goals(:,1),goals(:,2),goals(:,3));
xlim([-6 6]);
ylim([-7 7]);
zlim([0 5]);
xlabel('x');
ylabel('y');
zlabel('z');
axis vis3d;
camproj('perspective');
hold off;



