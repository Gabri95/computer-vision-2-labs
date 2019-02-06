

function [weights] = no_weights(source,target, source_normals, target_normals)
    [d, Ns] = size(source);
    weights = ones(1, Ns);
end

