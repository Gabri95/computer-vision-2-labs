function [res] = eight_point_algorithm(I1, I2, type, visualize, threshold)

    if isempty(visualize)
        visualize = false;
    end
    
    if isempty(threshold)
        threshold = 3;
    end

    % Compute matching between pictures using VLFeat
    [matches, feature1, feature2,desc1,desc2,score] = keypoint_matching(I1, I2, threshold);
    
    % Generate matrices p1 and p2 in format N_POINTS * [x, y, 1]
    p1 = [];
    p2 = [];
    for i = 1:1:size(matches,2)
       x1 = feature1(1, matches(1,i));
       x2 = feature2(1, matches(2,i));
       y1 = feature1(2, matches(1,i));
       y2 = feature2(2, matches(2,i));

       p1 = [p1; [x1, y1, 1]];
       p2 = [p2; [x2, y2, 1]];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    inliers_rateo = 0;
    time_to_best = 0;
    number_of_matches = 0;
    
    % Compute the best transformation F
    if type == "ransac" 
        [F, inliers, inliers_rateo, time_to_best, number_of_matches] = fundamental_matrix_RANSAC(p1, p2, false);
%         desc1 = desc1(:, inliers);
%         desc2 = desc2(:, inliers);
        desc1 = p1(inliers, :);
        desc2 = p2(inliers, :);
%         feature1 = feature1(:, inliers);
%         feature2 = feature2(:, inliers);
    elseif type == "normalized" 
        F = fundamental_matrix(p1, p2, true);
    else
        F = fundamental_matrix(p1, p2, false);
    end

    res.F = F;
    res.desc1 = desc1;
    res.desc2 = desc2;
    res.points1 = feature1;
    res.points2 = feature2;
    res.inliers_rateo = inliers_rateo;
    res.time_to_best = time_to_best;
    res.number_of_matches = number_of_matches;
    
    
    if visualize
        draw_epipolar_lines(p1, p2, I1, I2, F)
    end
    
end

