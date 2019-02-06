
function [target_matches, target_matches_normals] = match_kd(source, target, source_normals, target_normals) %, tree)
    [d, Nt] = size(target);
    
    tree = KDTreeSearcher(target');
    
    idxs = knnsearch(tree,source');
    target_matches = tree.X(idxs,:)';
    target_matches_normals = target_normals(:, idxs);
end


