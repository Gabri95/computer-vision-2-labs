close all
clear all
clc

addpath('./icp_variants/');

OUTPUT_FOLDER = './output0';

target = load(fullfile('Data','r_target.mat'));
source = load(fullfile('Data','r_source.mat'));
target = target.target;
source = source.source;


ss = load(fullfile('Data','r_source_normals.mat'));
source_normals = ss.source_normals;
ss = load(fullfile('Data','r_target_normals.mat'));
target_normals = ss.target_normals;

ss= load(fullfile('Data','rotation.mat'));
truth_angles = ss.rotate;

ss= load(fullfile('Data','translate.mat'));
truth_translate = ss.translation;

figure()
plotPCL(target, []);
hold on
plotPCL(source, []);
hold off
%return
source_normals = compute_normal_from_mesh(source);
target_normals = compute_normal_from_mesh(target);

samplings = {{'NoSampling', []},
             {'Normalspace', @(s, t, s_n, t_n) normalspace_sampling(s, t, s_n, t_n, 0.7, 0.7)},
             {'Uniform', @(s, t, s_n, t_n) uniform_sampling(s, t, s_n, t_n, 0.7, 0.7)}};

selections = {{'NoSelection', []},
             {'Random', @(s, t, s_n, t_n) random_selection(s, t, s_n, t_n, 0.7, 0.7)}};

matchings = {{'KDTree', @match_kd}};

rejections = {{'NoRejection', []},
              {'Dmax=0.5', 0.5},
              {'Dmax=1', 1},
              {'Dmax=2', 2},
              {'OrientationConsistency0.5', @(s, t, s_n, t_n) rejection_by_normals(s, t, s_n, t_n, 0.5)},
              {'OrientationConsistency0.8', @(s, t, s_n, t_n) rejection_by_normals(s, t, s_n, t_n, 0.8)},
              };
          
weightings = {{'NoWeights', @no_weights},
              {'Godin94', @weighting_Godin94},
              {'OrientationConsistencyW', @weighting_normals},
              }; 
sampling = samplings{1};
selection = selections{1};
matching = matchings{1};
rejection = rejections{1};
weighting = weightings{1};

for i=1:length(samplings)
    new_sampling = samplings{i};
    new_error = run_experiment(OUTPUT_FOLDER,source, target, source_normals, target_normals, new_sampling, selection, matching, rejection, weighting, false, truth_translate, truth_angles);
end

sampling = samplings{3};

for i=1:length(selections)
    new_selection = selections{i};    
    new_error = run_experiment(OUTPUT_FOLDER,source, target, source_normals, target_normals, sampling, new_selection, matching, rejection, weighting, false, truth_translate, truth_angles);
end


for i=1:length(rejections)
    new_rejection = rejections{i};    
    new_error = run_experiment(OUTPUT_FOLDER,source, target, source_normals, target_normals, sampling, selection, matching, new_rejection, weighting, false, truth_translate, truth_angles);
end


for i=3:length(weightings)
    new_weighting = weightings{i};    
    new_error = run_experiment(OUTPUT_FOLDER,source, target, source_normals, target_normals, sampling, selection, matching, rejection, new_weighting, false,truth_translate, truth_angles);
end
