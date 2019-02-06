
function [selected_sources, target_matches, source_normals, target_normals] = rejection_by_normals(selected_sources, target_matches, source_normals, target_normals, threshold)
    
    if isempty(threshold)
        threshold = 0.9;
    end
    
    dotproducts = sum(source_normals .* target_normals);
    
    idxs = dotproducts > threshold;
    
    selected_sources = selected_sources(:, idxs);
    target_matches = target_matches(:, idxs);
    source_normals = source_normals(:, idxs);
    target_normals = target_normals(:, idxs);
    
end
