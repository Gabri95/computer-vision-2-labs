clear all;
close all;

% Parse images

BASE_PATH = '../Data/House/';
names = dir(strcat(BASE_PATH, '*.png'));
images = cell(length(names), 1);
for i = 1:length(names)
   images{i} = imread(strcat(BASE_PATH, names(i).name));
end


% Run eight point algorithm
% type = 'basic';
% type = 'normalized';
type = 'ransac';

% thresholds = [0, 0.1, 0.5, 1, 2.5, 5, 7.5, 10, 12.5];
thresholds = [10];

for t=thresholds
    res = cell(1,length(names));
    for i = 2:size(images,1)
       res{i-1} = eight_point_algorithm(images{i-1}, images{i}, type, true, t);
    end

    total_time = 0;
    average_rateo = 0;
    average_matches = 0;

    for i = 1:48
       total_time = total_time + res{1,i}.time_to_best ;
       average_rateo = average_rateo + res{1,i}.inliers_rateo;
       average_matches = average_matches + res{1,i}.number_of_matches;
    end

    total_time = total_time/48;
    average_rateo = average_rateo/48;
    average_matches = average_matches/48;
    
    % Only meaningful if using RANSAC
    if strcmp(type,'ransac')
        disp({"Threshold", t, "Iterations", total_time, "AVG Rateo", average_rateo, "AVG_matches", average_matches}); 
    end
end
    