close all

living  =  load(fullfile('Data','livingRoom.mat'));
target = living.livingRoomData{1};
source = living.livingRoomData{3};
gridSize = 0.01;
fixed = pcdownsample(target, 'gridAverage', gridSize);
move = pcdownsample(source, 'gridAverage', gridSize);
% , 'MaxIterations',50,'Metric','pointToPlane'
T = pcregrigid(move,fixed,'Extrapolate',true);
T.T(4,1:3) = T.T(4,1:3) + [0.2,0,0];
move2 = pctransform(move,T);

target = double(fixed.Location');
target_normals = pcnormals(fixed)';

% target = pcread('teapot.ply');

A = eye(4,4);
Rot = eul2tform([-pi/10 pi/6 pi/8]);
%Rot = eul2tform([0 0 0]);
translation = [2,0,0.3];
%translation = [0,0,0];
A(4,1:3) = translation;
A(1:3,1:3) = Rot(1:3,1:3)

rotate = decompose_rotation(A(1:3,1:3));
tform1 = affine3d(A);
source_pcl = pctransform(move2,tform1)';

source = double(source_pcl.Location');
%source = source(1:3,30000:60000);
source_normals = pcnormals(source_pcl)';
%source_normals = source_normals(1:3,30000:60000);

save(fullfile('Data','r_target.mat'),'target');
save(fullfile('Data','r_target_normals.mat'),'target_normals');

save(fullfile('Data','r_source.mat'),'source');
save(fullfile('Data','r_source_normals.mat'),'source_normals');

save(fullfile('Data','rotation.mat'),'rotate');
save(fullfile('Data','translate.mat'),'translation');

figure()
plotPCL(target,[]);
hold on
plotPCL(source,[]);
hold off
