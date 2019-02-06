
function [weights] = weighting_normals(selected_sources, target_matches, source_normals, target_normals)
    
    
    weights = sum(source_normals .* target_normals);
end


