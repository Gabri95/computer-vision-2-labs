
function [source, target, source_normals, target_normals] = drop_borders(source, target, source_normals, target_normals, p_s, p_t)
   [source, source_normals] = drop_cloud_borders(source, source_normals, p_s);
   [target, target_normals] = drop_cloud_borders(target, target_normals, p_t);
end


function [points, normals] = drop_cloud_borders(points, normals, p)
    
    [d, N] = size(points);
    
    n = uint32(round(p*N));
    
    X_sorted = sort(squeeze(points(1, :)));
    Y_sorted = sort(squeeze(points(2, :)));
    
    x_min = X_sorted(1 + n);
    x_max = X_sorted(end - n);
    
    y_min = Y_sorted(1 + n);
    y_max = Y_sorted(end - n);
    
    idxs = (points(1, :) >= x_min) & (points(1, :) <= x_max) & (points(2, :) >= y_min) & (points(2, :) <= y_max);
    
    
    points = points(:, idxs);
    normals = normals(:, idxs);
    


end


