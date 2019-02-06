
function [selected_sources, target_matches, selected_source_normals, selected_target_normals] = random_selection(selected_sources, target_matches, source_normals, target_normals, p_s, p_t)
    
    [d, Nt] = size(target_matches);
    [d, Ns] = size(selected_sources);
    
    idxs = randi(Ns, 1, uint32(p_s*Ns));
    selected_sources = selected_sources(:, idxs);
    selected_source_normals = source_normals(:, idxs);
    
    idxs = randi(Nt, 1, uint32(p_t*Nt));
    target_matches = target_matches(:, idxs);
    selected_target_normals = target_normals(:, idxs);

    
    
%     idxs = (rand(1, Ns) - p_s) <= 0;
%     selected_sources = selected_sources(:, idxs);
%     selected_source_normals = source_normals(:, idxs);
%     
%     idxs = (rand(1, Nt) - p_t) <= 0;
%     target_matches = target_matches(:, idxs);
%     selected_target_normals = target_normals(:, idxs);
    
end

