function [POINTS] = structure_from_motion(PVM, subset_size, build_blocks, eliminate_affinity, average)

    [m, n] = size(PVM);
    
    assert (mod(m, 2) == 0);
    
    m = m / 2;
    
    connections = compute_connections(PVM);
    
    if isempty(build_blocks)
        build_blocks = build_chain;
    end
    
    Blocks = build_blocks(connections, subset_size);
   
    
    POINTS = zeros(3, n);
    points_count = zeros(1, n);
    
    points_set = false(1, n);
    
    for i = 1:size(Blocks, 1)
        block = Blocks(i, :);
        
        views = block{1};
        
        points = block{2};

        
        used_points = sum(points);
        
        
        rows = sort([2*views-1, 2*views]);
        
        
        assert(~any(any(isnan(PVM(rows, points)))));
        
        old_points_idxs = points_set(points);
        new_points_idxs = ~old_points_idxs;
        
%         figure(14);
%         pc = reshape(PVM(rows, points), 2, []);
%         colors = 1:size(pc, 2);
%         colors = 1 + uint8(colors/length(views));
%         scatter(pc(1, :)', pc(2, :)', 4, colors, 'filled')
%         drawnow
%         pause(1)
        
        
        if average || any(new_points_idxs)
        
            [M, S] = factorization(PVM(rows, points), eliminate_affinity);
            
            f = figure(14);
            scatter3(S(1, :)', S(2, :)', S(3, :)', 4, 'filled');
            title(sprintf('Reconstruction From Block %d', i));
            movegui(f, 'north');
            drawnow
%             pause(2)
            
            if ~any(points_set)
                POINTS(:, points) = S;
                plot_dense_block(connections, views, points, points_set, m, n)
                fprintf("%d points used;\n\n", used_points);
            else
    
                fprintf("%d points used; %d new points\n\n", used_points, sum(new_points_idxs));

                plot_dense_block(connections, views, points, points_set, m, n)

                [d,Z,transform] = procrustes(POINTS(:, points_set & points)', S(:, old_points_idxs)'); %, 'reflection', false);
                Z = transform.b*S'*transform.T + transform.c(1, :);
                Z = Z';

%                 fprintf("(%d, %d), (%d, %d), (%d, %d)\n", size(transform.b, 1), size(transform.b, 2), size(transform.T, 1), size(transform.T, 2), size(transform.c, 1), size(transform.c, 2));
%                 fprintf("%f, (%f, %f, %f)\n", transform.b, transform.c(1, 1), transform.c(1, 2), transform.c(1, 3));

                if average
                    %%% take a weighted average of the estimations, with weights
                    %%% proportional to the number of camera used to evaluate them

                    POINTS(:, points) = POINTS(:, points ).*points_count(points) +  length(views)*Z;
                    points_count = points_count + length(views)*points;
                    POINTS(:, points) = POINTS(:, points )./points_count(points);
                else
                    %%% give preference to estimations coming from blocks closer
                    %%% to the root
                    POINTS(:, points & ~points_set) = Z(:, new_points_idxs);
                end
                
                f=figure(15);
                
                scatter3(S(1, old_points_idxs)', S(2, old_points_idxs)', S(3, old_points_idxs)', 8, 'red', 'filled')
                hold on
                scatter3(POINTS(1, points_set & points)', POINTS(2, points_set & points)', POINTS(3, points_set & points)', 8, 'green', 'filled')
                scatter3(Z(1, old_points_idxs)', Z(2, old_points_idxs)', Z(3, old_points_idxs)', 8, 'blue', 'filled')
                hold off
                
                legend(["Recontruction", "Base Point Cloud", "Transformed"]);
                
                title(sprintf('Transformation For Block %d', i));
                
                movegui(f, 'south');
                drawnow
            end
            
            colors = zeros(length(points_set), 3);
            colors(points & points_set, 3) = 1;
            colors(points & ~points_set, 1) = 1;

            points_set = points_set | points;

            f = figure(10);
            scatter3(POINTS(1, points_set)', POINTS(2, points_set)', POINTS(3, points_set)', 4, colors(points_set, :), 'filled');
            title('Current Reconstruction');
            movegui(f, 'east');
            drawnow

            
        end
        
        
    end
    
    POINTS = POINTS(:, points_set);
    
end


function plot_dense_block(connections, views, points, points_set, m, n)
        bin_mat = 0.8*ones(m, n, 3);
%         bin_mat(:, points_set, [2, 3]) = 0;
        bin_mat([connections, connections, connections]) = 0;
        bin_mat(views, points & ~points_set, 3) = 1;
        bin_mat(views, points & points_set, 1) = 1;
        
        f=figure(12);
        imshow(bin_mat,'InitialMagnification','fit');
        daspect([n,m,1]);
        movegui(f, 'west');
        drawnow
end



