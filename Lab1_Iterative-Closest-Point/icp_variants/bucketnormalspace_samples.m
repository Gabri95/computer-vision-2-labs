
function [matches, matches_normals] = bucketnormalspace_samples(points, normals, p)

    [d, N] = size(points);
    
    if isfloat(p)
        p = round(p*N);
    end
    
    
    BINS = 648;
    bins_centers = generate_sphere_bins(BINS);
%     plotPCL(bins_centers, []);

    bin_tree = KDTreeSearcher(bins_centers');

    matchings = knnsearch(bin_tree, normals');
    
    nonempty_bins = unique(matchings);
    
    matchings = [matchings, [1:size(normals, 2)]'];
    
    
    BINS = size(nonempty_bins, 1);
    
    sampled_buckets = mnrnd(p, ones(1, BINS)/BINS, 1); 
    
    matches_normals = [];
    matches = [];
    
    
    for bin=1:BINS
        b = nonempty_bins(bin);
        idxs = matchings(:, 1) == b;
        matched = matchings(idxs, 2);
        
        
        dim = size(matched, 1);
        
        if dim > 0
            normals_subset = normals(:, matched);
            points_subset = points(:, matched);

            samples = randi(dim, 1, sampled_buckets(bin));

            matches = [matches, points_subset(:, samples)];
            matches_normals = [matches_normals, normals_subset(:, samples)];
        else
            disp('I never appear');
        end
    end

    
%     matchings = bins_centers' * normals;
%     
%     matchings = matchings == max(matchings); 
%     
%     nonempty_bins = sum(matchings, 2) > 0;
%     matchings = matchings(nonempty_bins, :);
%     
%     BINS = size(matchings, 1);
%     
%     sampled_buckets = mnrnd(p, ones(1, BINS)/BINS, 1); 
%     
%     matches_normals = [];
%     matches = [];
%     
%     
%     for b=1:BINS
%         dim = sum(matchings(b, :));
%         if dim > 0
%             normals_subset = normals(:, matchings(b, :));
%             points_subset = points(:, matchings(b, :));
% 
%             samples = randi(dim, 1, sampled_buckets(b));
% 
%             matches = [matches, points_subset(:, samples)];
%             matches_normals = [matches_normals, normals_subset(:, samples)];
%         end
%     end
    
    
    
end



