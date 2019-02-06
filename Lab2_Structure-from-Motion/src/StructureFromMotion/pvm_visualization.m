clc
clear all;
close all;


BASE_PATH = '../Data/House/';

% Parse images

names = dir(strcat(BASE_PATH, '*.png'));

N_VIEWS = length(names);

images = cell(N_VIEWS, 1);
for i = 1:length(names)
   images{i} = imread(strcat(BASE_PATH, names(i).name));
end


step = 8;

use_supplied_data = false;



if use_supplied_data
    
    D = load('../PointViewMatrix.txt');
    [m, n] = size(D);
    
%     D = D(1:(m-mod(m, 4))/2, :); %1:uint16(n/2));
    [m, n] = size(D);
    m = m / 2;
    r = uint16(n/m);
    conn = false(m,n);

    span = step;
    for i =1:m
        conn(i:min(i, end), r*(i-1)+1:min(r*(i + step+ randi(span)), end)) = true;
    end
    conn(end, r*i:end) = true;
    
    rows = 1:m;
    
    for j =1:n
        i = max(rows(conn(:, j)));
        d = uint16(0.7*(m - i));
        if d > 0 && randi(3) > 1
            conn(i:i+randi(d), j) = true;
        end
    end
%     
%     for i =1:m
%         D(2*i-1:2*i, ~conn(i, :)') = NaN;
%     end
    
else
    
    D = load('../pointview.mat', '-ASCII');
%     D(isnan(D)) = 3000;
%     orig_size = size(D);
%     D = reshape(D, [], 1);
%     vals = 100*randn(sum(isnan(D)), 1);
%     
%     D(isnan(D)) = vals; 
%     D = reshape(D, orig_size);
end

conn = compute_connections(D);

pause(1)

[m, n] = size(D);


for i=3:2:m
    figure(2);
    idx = (i+1)/2;
    imshow(images{idx});
    title(idx);
    hold on
%     scatter(D(i-2, ~isnan(D(i, :)))', D(i-1, ~isnan(D(i, :)))', 8, 'blue', 'filled');
%     scatter(D(i, :)', D(i+1, :)', 8, 'red', 'filled');
%     scatter(D(i-2, isnan(D(i, :)))', D(i-1, isnan(D(i, :)))', 8, 'yellow','filled');
    
    scatter(D(i-2, :)', D(i-1, :)', 8, 'red', 'filled');
    scatter(D(i, :)', D(i+1, :)', 8, 'yellow', 'filled');
    
    matches = D(i-2:i+1, :);
    matches = matches(:, ~any(isnan(matches)));
    plot(matches([1, 3], :), matches([2, 4], :));
% 
% 
    hold off
    axis equal
    pause(1)
end
