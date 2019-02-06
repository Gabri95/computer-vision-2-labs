function [Blocks] = find_maximal_cliques_04(connections, Nv, Np)

% from Sec. 3.2.2 of the paper "ijcv04d.pdf" provided in "Reading"

[m, n] = size(connections);

views = 1:m;

Blocks =  [];

for i =1:n
    
    rows = views(connections(:, i));
    points = all(connections(rows, :));

    if length(rows) > Nv && sum(points) > Np
        Block = {rows, points};
        Blocks = [Blocks; Block];
    end
    
end

end
