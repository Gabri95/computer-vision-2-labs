function plotPCL(A, color)
    X = A(1,:);
    Y = A(2,:);
    Z = A(3,:);
    axis equal
    if isempty(color)
        scatter3(X,Y,Z, 1, 'filled');
    else
        scatter3(X,Y,Z, 1, color, 'filled');
    end
    axis equal
end