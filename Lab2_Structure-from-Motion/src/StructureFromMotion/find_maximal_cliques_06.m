function [Blocks] = find_maximal_cliques_06(connections, Nv, Np)

%from APPENDIX of http://www-cvr.ai.uiuc.edu/ponce_grp/publication/paper/pami06.pdf

[m, n] = size(connections);

last = ones(1, n);
first = m * ones(1, n);

%here we assume that each point appears in a consecutive set of views in
%the PVM matrix. Therefore, a column will have the form:
%[NaN, NaN, ... , NaN, x1, x2, ..., xm, NaN, NaN, ..., NaN]


for i =1:n
   j = 1;
   while ~connections(j, i) && j < m
       j = j + 1;
   end
   while connections(j, i) && j < m
       j = j + 1;
   end
   while j < m
       assert(~connections(j, i));
       j = j + 1;
   end
end

for i =1:n
   j = 1;
   while ~connections(j, i) && j < m
       j = j + 1;
   end
   while connections(j, i) && j < m
       j = j + 1;
   end
   while j < m
       assert(~connections(j, i));
       j = j + 1;
   end
end


views = 1:m;
for i =1:n
    point_views = views(connections(:, i));
    first(i) = min(point_views);
    last(i) = max(point_views);
    
end

tracks = 1:n;

tracks = tracks(last - first - (Nv - 1)> 0);

starting_views = sort(unique(first(tracks)));
ending_views = sort(unique(last(tracks)));


Blocks =  [];

for ii =1:length(starting_views)
    i = starting_views(ii);
    for jj =1:length(ending_views)
        j = ending_views(jj);
        
        if j > i
            B = connections(i, :) & connections(j, :);
            
            if ismember(i, first(B)) && ismember(j, last(B)) && sum(B) > Np

                assert (all(all(connections(i:j, B))));
                
                Block = {i:j, B};

                Blocks = [Blocks; Block];
            end
            
        end
    end
end

end
