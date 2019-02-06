
function [Blocks] = build_tree(connections, subset_size, build_blocks)
    % Ispired from
    % http://www-cvr.ai.uiuc.edu/ponce_grp/publication/paper/pami06.pdf (for the initial dense block) and
    % "3d object modeling and recognition using local affine-invariant image descriptors and multi-view spatial constraints"
    % for the "stitch graph", though here we use a Maximum Spanning Tree instead of Dijkstra Algorithm
    
    
    if isempty(build_blocks)
        build_blocks = build_chain;
    end
    
    Blocks = build_blocks(connections, subset_size);
    
%     Blocks = find_maximal_cliques_06(connections, subset_size, 5);
%     Blocks = find_maximal_cliques_04(connections, subset_size, 5);
%     Blocks = build_chain(connections, subset_size);

    if size(Blocks, 1) > 1
        
        score = -1;
        root = -1;

        for i= 1:size(Blocks, 1)
           if score < sum(Blocks{i, 2})%*length(Blocks{i, 1})
               score = sum(Blocks{i, 2});%*length(Blocks{i, 1});
               root = i;
           end
        end

        s = [];
        t = [];
        w = [];

        for i= 1:size(Blocks, 1)-1
           for j= i+1:size(Blocks, 1)
               if sum(Blocks{i, 2} & Blocks{j, 2}) > 0
                   s = [s, i];
                   t = [t, j];
                   w = [w, -1*sum(Blocks{i, 2} & Blocks{j, 2})];
               end
           end
        end

        G = graph(s,t,w);

        [Tree,pred] = minspantree(G, 'Root', root);


        connected = ~isnan(pred) & (pred ~= 0);
        nodes = 1:Tree.numnodes;

        DirectedTree = digraph(pred(connected), nodes(connected));

        DirectedTree = rmnode(DirectedTree, nodes(isnan(pred)));
% 
%         figure;
%         plot(DirectedTree);


        order = toposort(DirectedTree);

        Blocks = Blocks(order, :);
    end
    
end

