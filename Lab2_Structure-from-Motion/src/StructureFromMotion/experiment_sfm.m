clc
clear all;
close all;


step = 3;

use_supplied_data = false;
sparsify_data = false;

remove_affinity = true;
average = true;

build_blocks = @(connections, subset_size) build_chain(connections, subset_size);
% build_blocks = @(connections, subset_size) build_tree(connections, subset_size, @build_chain);
% build_blocks = @(connections, subset_size) build_tree(connections, subset_size, @(c, s) find_maximal_cliques_06(c, s, 5));
% build_blocks = @(connections, subset_size) build_tree(connections, subset_size, @(c, s) find_maximal_cliques_04(c, s, 10));


if use_supplied_data
    
    D = load('../PointViewMatrix.txt');
    [m, n] = size(D);
    
    
    if sparsify_data
    %     D = D(1:(m-mod(m, 4))/2, :); %1:uint16(n/2));
        [m, n] = size(D);
        m = m / 2;
        r = uint16(n/m);
        conn = false(m,n);

        span = step;
        for i =1:m
            conn(i:min(i, end), r*(i-1)+1:min(r*(i+ step +randi(span)), end)) = true;
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
        
        for i =1:m
            D(2*i-1:2*i, ~conn(i, :)') = NaN;
        end
    end
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




if use_supplied_data && ~sparsify_data
    [M, S] = factorization(D, remove_affinity);
else
    S = structure_from_motion(D, step, build_blocks, remove_affinity, average);
end

figure;
scatter3(S(1, :)', S(2, :)', -S(3, :)', 4, 'filled');
% axis equal

