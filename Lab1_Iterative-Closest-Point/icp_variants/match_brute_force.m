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

