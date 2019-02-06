
function [matches, matches_normals] = uniformnormalspace_samples(points, normals, p)

    [d, N] = size(points);
    
    if isfloat(p)
        p = uint32(p*N);
    end
    
    tree = KDTreeSearcher(normals');
    
    
    samples = randn(3, p);
    samples = samples ./ sum(samples.^2);
    
    
    idxs = knnsearch(tree, samples');
    
    matches_normals = normals(:, idxs);
    matches = points(:, idxs);
end

