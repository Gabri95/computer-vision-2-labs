function [transformed_pts, Rotation, Translation, final_error, iterations] = icp(source, target, source_normals, target_normals, sampling, selection, matching, rejection, weighting, show_plot)
    
    tic;
    
    if isempty(show_plot)
        show_plot = false;
    end
    
    if isempty(sampling)
        sampling = @no_selection;
    end
    
    if isempty(selection)
        selection = @no_selection;
    end
    
    if isempty(weighting)
        weighting = @no_weights;
    end
    
    D_max = -1;
    
    if isnumeric(rejection)
        D_max = rejection;
    elseif isempty(rejection)
        D_max = -1;
    elseif ~ isa(rejection, 'function_handle')
        disp('ERROR')
    end
    
    
    [source, target, source_normals, target_normals] = sampling(source,target,source_normals,target_normals);
    
    [d, Nt] = size(target);
    
    source_pts = source;
    
    
    dist = 0;
    % error_trend = [rms(source_pts, target)];
    [matched_target, ~] = matching(source_pts, target, source_normals, target_normals);
    old_dist = rms(source_pts, matched_target);
    error_trend = [old_dist];
    
    if show_plot
        fig = figure;
        error_trend_plot = plot(error_trend);
        xlabel('Iteration')
        ylabel('Distance')
    end
    
    Rotation = eye(3);
    Translation = [0; 0; 0];
    iters = 1;
    while abs(old_dist - dist)/(old_dist) > 0.0001 && iters < 1000 %abs(dist - old_dist) > 0.000001
        iters = iters + 1;
        [selected_source, selected_target, selected_source_normals, selected_target_normals] = selection(source_pts,target, source_normals, target_normals);
        
        [target_matches, target_matches_normals] = matching(selected_source, selected_target, selected_source_normals, selected_target_normals);

        if D_max >= 0
            [selected_source, target_matches, selected_source_normals, target_matches_normals] = rejection_by_distance(selected_source, target_matches, selected_source_normals, target_matches_normals, D_max);
        elseif isa(rejection, 'function_handle')
            [selected_source, target_matches, selected_source_normals, target_matches_normals] = rejection(selected_source, target_matches, selected_source_normals, target_matches_normals);
        end
        
        weights = weighting(selected_source, target_matches, selected_source_normals, target_matches_normals);
        
        if D_max >= 0
            distances = sum((selected_source - target_matches).^2);
            D_max = update_threshold(distances, D_max);
        end
        
        target_centroid = centroid(target_matches,weights);
        source_centroid = centroid(selected_source,weights);
        
        Y = target_matches - target_centroid;
        X = selected_source - source_centroid;
        
        %Building the diagonal matrix is too space-expensive (quadratic in
        %the number of points) and can't work for "merge scene" as we can 
        %have more that 50.000 points
        
        %S = X * diag(weights) * Y';
        S = X* (Y'.*weights');
        
        [U,~,V] = svd(S);
        M = eye(d, d);
        M(d, d) = det(V*U');
        
        R = V * M * U';
        trans = target_centroid - R*source_centroid;
        
        source_pts = R*source_pts + trans;
        source_normals = R*source_normals;
        
        Rotation = R*Rotation;
        Translation = R*Translation + trans;
        
        old_dist = dist;
        
        [matched_target, ~] = matching(source_pts, target, source_normals, target_normals);
        dist = rms(source_pts, matched_target);
        
        error_trend = [error_trend, dist];
        
        if show_plot
            set(error_trend_plot,'XData',1:length(error_trend))
            set(error_trend_plot,'YData',error_trend)
            refreshdata
            drawnow
        end
    end
    
    fprintf('%f seconds elapsed\n', toc);
    fprintf('Final distance: %f\n', dist);
    fprintf('Number of Iterations: %d\n', length(error_trend));
    
    transformed_pts = source_pts;
    
    iterations = length(error_trend);
    final_error = error_trend(end);
    
end


function [score] = rms(source,target)
    score = sum(norm(source - target));
end

function [center] = centroid(points, weights)
    weighted_center = sum(points .* weights,2);
    center =  weighted_center / sum(weights);
end

function [center] = centroid_median(points, weights)
%     tot_weight = sum(weights);
    
    center = [median(points(1, :)); median(points(2, :)); median(points(3, :))];
    
end


function [D_max] = update_threshold(distances, D)
    
    mu = mean(distances);
    sigma = std(distances);
    
    if mu < D
        D_max = mu + 3*sigma;
    elseif mu < 3*D
        D_max = mu + 2*sigma;
    elseif mu < 6*D
        D_max = mu + sigma;
    else
        if(isempty(distances))
            disp("no points")
            D_max = 1;
        else
            [n, d] = hist(distances);
        
            [highest_peak, bin] = max(n); 
        
            D_max = d(bin);
            bin = bin + 1;
        
            while bin <= lenght(n) && n(bin) > 0.6*highest_peak
                D_max = d(bin);
                bin = bin + 1;
            end
        end
    end
end


function [selected_sources, target_matches, source_normals, target_normals] = rejection_by_distance(selected_sources, target_matches, source_normals, target_normals, D_max)
    
    distances = sum((selected_sources - target_matches).^2);
    
    idxs = distances < D_max;
    
    selected_sources = selected_sources(:, idxs);
    target_matches = target_matches(:, idxs);
    source_normals = source_normals(:, idxs);
    target_normals = target_normals(:, idxs);
    
end

function [sources, targets, source_normals, target_normals] = no_selection(sources, targets, source_normals, target_normals)
end

function [weights] = no_weights(sources,targets, source_normals, target_normals)
    [d, Ns] = size(sources);
    weights = ones(1, Ns);
end
