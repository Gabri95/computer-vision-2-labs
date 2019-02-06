
function [selected_sources, target_matches, source_normals, target_normals] = rejection_by_normals(selected_sources, target_matches, source_normals, target_normals)
    
    dotproducts = sum(source_normals .* target_normals);
    
    idxs = dotproducts > 0.95;
    
    selected_sources = selected_sources(:, idxs);
    target_matches = target_matches(:, idxs);
    source_normals = source_normals(:, idxs);
    target_normals = target_normals(:, idxs);
    
end



function [weights] = no_weights(source,target, source_normals, target_normals)
    [d, Ns] = size(source);
    weights = ones(1, Ns);
end

function [weights] = weighting_Godin94(selected_sources, target_matches, source_normals, target_normals)
    distances = sum((selected_sources - target_matches).^2);
    
    max_distance = max(distances);
    
    weights = 1 - (distances / max_distance);
end

function [weights] = weighting_normals(selected_sources, target_matches, source_normals, target_normals)
    
    
    weights = sum(source_normals .* target_normals);
end

function selected_sources = match_brute_force(source,target)
    [d, Ns] = size(source);
    [d, Nt] = size(target);
    
    target_matches = zeros(d,Ns);
    
    for s=1:Ns
        target_matches(:,s) = target(:,1);
        min_dist = 10000000;
        for t=1:Nt
            res = norm(target(:,t) - source(:,s));
            if  res < min_dist
                min_dist = res;
                target_matches(:,s) = target(:,t);
            end
        end
    end
end

function [target_matches, target_matches_normals] = match_kd(source, target, source_normals, target_normals) %, tree)
    [d, Nt] = size(target);
    
    tree = KDTreeSearcher(target');
    
    idxs = knnsearch(tree,source');
    target_matches = tree.X(idxs,:)';
    target_matches_normals = target_normals(:, idxs);
end


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


function [points, normals] = drop_borders(points, normals, p)
    
    [d, N] = size(points);
    
    n = uint32(round(p*N));
    
    X_sorted = sort(squeeze(points(1, :)));
    Y_sorted = sort(squeeze(points(2, :)));
    
    x_min = X_sorted(1 + n);
    x_max = X_sorted(end - n);
    
    y_min = Y_sorted(1 + n);
    y_max = Y_sorted(end - n);
    
    idxs = (points(1, :) >= x_min) & (points(1, :) <= x_max) & (points(2, :) >= y_min) & (points(2, :) <= y_max);
    
    
    points = points(:, idxs);
    normals = normals(:, idxs);
    


end


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


function vectors = generate_sphere_bins(N)
%     https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
%     vectors = zeros(3, N);
% 
%     count = 0;
%     a = 4 * pi/N;
%     d = sqrt(a);
% 
%     M_theta = round(pi/d);
%     d_theta = pi/M_theta;
%     d_phi = a/ d_theta;
%     for m = 0:M_theta-1
%         theta = pi*(m + 0.5)/M_theta;
%         M_phi = round(2*pi*sin(theta/d_phi));
%         for n =0:M_phi - 1
%             phi = 2*pi*n / M_phi;
%             vectors(:, count + 1) = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
%             count = count + 1;
%         end
%     end
%     count


    vectors = zeros(3, N);
    offset = 2./N;
    increment = pi * (3 - sqrt(5));
    rnd = 1; %randi(N);
    
    for i = 0:N-1
        y = ((i * offset) - 1) + (offset / 2);
        r = sqrt(1 - y^2);

        phi = mod(i + rnd, N) * increment;

        x = cos(phi) * r;
        z = sin(phi) * r;

        vectors(:, i+1) = [x;y;z];
    end
end

function [matches, matches_normals] = bucketnormalspace_samples(points, normals, p)

    [d, N] = size(points);
    
    if isfloat(p)
        p = round(p*N);
    end
    
    
    BINS = 648;
    bins_centers = generate_sphere_bins(BINS);
%     plotPCL(bins_centers, []);
    matchings = bins_centers' * normals;
    
    matchings = matchings == max(matchings); 
    
    nonempty_bins = sum(matchings, 2) > 0;
    matchings = matchings(nonempty_bins, :);
    
    BINS = size(matchings, 1);
    
    sampled_buckets = mnrnd(p, ones(1, BINS)/BINS, 1); 
    
    matches_normals = [];
    matches = [];
    
    
    for b=1:BINS
        dim = sum(matchings(b, :));
        if dim > 0
            normals_subset = normals(:, matchings(b, :));
            points_subset = points(:, matchings(b, :));

            samples = randi(dim, 1, sampled_buckets(b));

            matches = [matches, points_subset(:, samples)];
            matches_normals = [matches_normals, normals_subset(:, samples)];
        end
    end
    
    
    
end


function [source, target, source_normals, target_normals] = normalspace_sampling(source, target, source_normals, target_normals, ps, pt)
   [source, source_normals] = bucketnormalspace_samples(source, source_normals, ps);
   [target, target_normals] = bucketnormalspace_samples(target, target_normals, pt);
end

