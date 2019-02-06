
function [Blocks] = build_chain(connections, subset_size)

    IMAGES = size(connections, 1);
    POINTS = size(connections, 2);
    Blocks = [];
    for i=1:(IMAGES - subset_size)
        mask = true(1, POINTS);
        for j=i:i+subset_size
            mask = mask & connections(j,:); 
        end
        Blocks = [Blocks; {i: i+subset_size, mask}];
    end

end

