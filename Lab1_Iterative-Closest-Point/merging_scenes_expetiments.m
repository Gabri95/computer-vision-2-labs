close all
clear all
clc

addpath('./icp_variants/');

%DATASET = './data_mat1/';
DATASET = './Data/data_mat1/';


files = dir(fullfile(DATASET, '*.mat'));
nfiles = length(files);

n_useful_files = 0;

for i = 1:nfiles
    [pathstr, name, ext] = fileparts(files(i).name);
    parts = strsplit(name, '_');
    
    if length(parts) == 1 && ext == "mat"
        n_useful_files = n_useful_files +1;
    end
end



frames = cell(n_useful_files); 
frames_normals = cell(n_useful_files);

for i = 1:nfiles
    
    [pathstr, name, ext] = fileparts(files(i).name);
    parts = strsplit(name, '_');
    
    if ext == ".mat"
        if length(parts) == 1
            
            frame = load(fullfile(DATASET, files(i).name));

            point_cloud = frame.points';

%             point_cloud =  point_cloud(:, point_cloud(3, :) < 1.5);

            frames{frame.id + 1} = point_cloud;
        elseif parts(2) == "normal"
            
            normals = load(fullfile(DATASET, files(i).name));

            normal = normals.normal';

            frames_normals{normals.id + 1} = normal;
            
            
        end
    end
end

for f=1:length(frames)
    frame = frames{f};
    normals = frames_normals{f};
    
    mask = frame(3, :) < 1.25;
    
    frames{f} = frame(:, mask);
    frames_normals{f} = normals(:, mask);
    
end



samplings = {{'NoSampling', []},
             {'Normalspace0.8', @(s, t, s_n, t_n) normalspace_sampling(s, t, s_n, t_n, 0.8, 0.8)},
             {'Normalspace0.5', @(s, t, s_n, t_n) normalspace_sampling(s, t, s_n, t_n, 0.5, 0.5)},
             {'Uniform0.5', @(s, t, s_n, t_n) uniform_sampling(s, t, s_n, t_n, 0.5, 0.5)},
             {'Uniform0.8', @(s, t, s_n, t_n) uniform_sampling(s, t, s_n, t_n, 0.8, 0.8)}};

selections = {{'NoSelection', []},
             {'Random0.3', @(s, t, s_n, t_n) random_selection(s, t, s_n, t_n, 0.3, 0.3)},
             {'Random0.5', @(s, t, s_n, t_n) random_selection(s, t, s_n, t_n, 0.5, 0.5)},
             {'Random0.7', @(s, t, s_n, t_n) random_selection(s, t, s_n, t_n, 0.7, 0.7)}};
         

matchings = {{'KDTree', @match_kd}};

rejections = {{'NoRejection', []},
              {'Dmax=0.5', 0.5},
              {'Dmax=1', 1},
              {'Dmax=2', 2},
              {'OrientationConsistency0.5', @(s, t, s_n, t_n) rejection_by_normals(s, t, s_n, t_n, 0.5)},
              {'OrientationConsistency0.75', @(s, t, s_n, t_n) rejection_by_normals(s, t, s_n, t_n, 0.75)},
              {'OrientationConsistency0.9', @(s, t, s_n, t_n) rejection_by_normals(s, t, s_n, t_n, 0.9)},
              };
          
weightings = {{'NoWeights', @no_weights},
              {'Godin94', @weighting_Godin94},
              {'OrientationConsistencyW', @weighting_normals},
              };          

rates = {5, 1, 2, 4, 10};

references = {{'FirstFrame', 1},
             {'MiddleFrame', []}};
         
merges = {false, true};


sampling = samplings{1};
selection = selections{1};
matching = matchings{1};
rejection = rejections{1};
weighting = weightings{1};
rate = rates{2};
reference = references{1};
merge = merges{1};

best_error = run_experiment(frames, frames_normals, rate, merge, reference, sampling, selection, matching, rejection, weighting, true);

% best_error = 10;
% 
% for i=5:length(rates)
%     new_rate = rates{i};
%     
%     new_error = run_experiment(frames, frames_normals, new_rate, merge, reference, sampling, selection, matching, rejection, weighting, true);
%     
%     if new_error < best_error
%         best_error = new_error;
%         rate = new_rate;
%     end 
% end
% 
% return

for i=2:length(samplings)
    new_sampling = samplings{i};
    new_error = run_experiment(frames, frames_normals, rate, merge, reference, new_sampling, selection, matching, rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        sampling = new_sampling;
    end 
end

for i=2:length(selections)
    new_selection = selections{i};
    
    new_error = run_experiment(frames, frames_normals, rate, merge, reference, sampling, new_selection, matching, rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        selection = new_selection;
    end 
end

for i=2:length(matchings)
    new_matching = matchings{i};
    
    new_error = run_experiment(frames, frames_normals, rate, merge, reference, sampling, selection, new_matching, rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        matching = new_matching;
    end 
end

for i=2:length(rejections)
    new_rejection = rejections{i};
    
    new_error = run_experiment(frames, frames_normals, rate, merge, reference, sampling, selection, matching, new_rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        rejection = new_rejection;
    end 
end

for i=2:length(weightings)
    new_weighting = weightings{i};
    
    new_error = run_experiment(frames, frames_normals, rate, merge, reference, sampling, selection, matching, rejection, new_weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        weighting = new_weighting;
    end 
end


for i=2:length(rates)
    new_rate = rates{i};
    
    new_error = run_experiment(frames, frames_normals, new_rate, merge, reference, sampling, selection, matching, rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        rate = new_rate;
    end 
end

for i=2:length(references)
    new_reference = references{i};
    
    new_error = run_experiment(frames, frames_normals, rate, merge, new_reference, sampling, selection, matching, rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        reference = new_reference;
    end 
end

for i=2:length(merges)
    new_merge = merges{i};
    
    new_error = run_experiment(frames, frames_normals, rate, new_merge, reference, sampling, selection, matching, rejection, weighting, true);
    
    if new_error < best_error
        best_error = new_error;
        merge = new_merge;
    end 
end





function [error] = run_experiment(frames, frames_normals, rate, iterative_merge, reference, sampling, selection, matching, rejection, weighting, plot)
    
    fprintf('\n\nRate: %d; Merging: %s; Reference: %s, Sampling: %s; Selection: %s; Matching: %s; Rejection: %s; Weighting: %s\n\n', rate, mat2str(iterative_merge), reference{1}, sampling{1}, selection{1}, matching{1}, rejection{1}, weighting{1});
    
    name = sprintf('./output2/%d_%s_%s_%s_%s_%s_%s_%s', rate, mat2str(iterative_merge), reference{1}, sampling{1}, selection{1}, matching{1}, rejection{1}, weighting{1});
    
    [cloud, cloud_normals, error] = merge_scene(frames, frames_normals, rate, iterative_merge, reference{2}, sampling{2}, selection{2}, matching{2}, rejection{2}, weighting{2}, plot, name);

    close all
end

