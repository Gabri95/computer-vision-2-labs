
function [source, target, source_normals, target_normals] = uniform_sampling(source, target, source_normals, target_normals, p_s, p_t)


    [d, Nt] = size(target);
    [d, Ns] = size(source);
    
    idxs = (rand(1, Ns) - p_s) <= 0;
    source = source(:, idxs);
    source_normals = source_normals(:, idxs);
    
    idxs = (rand(1, Nt) - p_t) <= 0;
    target = target(:, idxs);
    target_normals = target_normals(:, idxs);
%       [source, target, source_normals, target_normals] = random_selection(source, target, source_normals, target_normals, p_s, p_t);    
end


