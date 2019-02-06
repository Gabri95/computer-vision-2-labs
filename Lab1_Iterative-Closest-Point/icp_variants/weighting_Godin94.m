
function [weights] = weighting_Godin94(selected_sources, target_matches, source_normals, target_normals)
    distances = sum((selected_sources - target_matches).^2);
    
    max_distance = max(distances);
    
    weights = 1 - (distances / max_distance);
end

