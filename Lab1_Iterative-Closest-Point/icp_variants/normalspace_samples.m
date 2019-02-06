
function [matches, matches_normals] = normalspace_samples(points, normals, p)

    [d, N] = size(points);
    
    if isfloat(p)
        p = uint32(p*N);
    end
    
    tree = KDTreeSearcher(normals');
    
    
    yaw = rand(1, p)* 2*pi;
    pitch = rand(1, p) * pi - 0.5*pi;
    
    samples = zeros(3, p);
    samples(1, :) = cos(yaw).*cos(pitch);
    samples(2, :) = sin(yaw).*cos(pitch);
    samples(3, :) = sin(pitch);
    
    idxs = knnsearch(tree, samples');
    
    matches_normals = normals(:, idxs);
    matches = points(:, idxs);
    
    
end


