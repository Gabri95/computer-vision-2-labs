
function [source, target, source_normals, target_normals] = normalspace_sampling(source, target, source_normals, target_normals, ps, pt)
   [source, source_normals] = bucketnormalspace_samples(source, source_normals, ps);
   [target, target_normals] = bucketnormalspace_samples(target, target_normals, pt);
end
