function draw_epipolar_lines(p1, p2, I1, I2, F)
    probability = 1;
    idxs = (rand(1, size(p2, 1)) - probability) <= 0;
    
    s = p1(idxs, :);
    S = I1;
    
    t = p2(idxs, :);
    T = I2;
    
    
    figure(1);
    
    subplot(1, 2, 1);
    imshow(S)
    hold on
    scatter(s(:,1), s(:, 2), 30, 'green', 'filled')
    hold off
    
    subplot(1, 2, 2);
    imshow(T)
    hold on
    for i = 1:size(s, 1)
        p = F*s(i,:)';
        
        x = linspace(0, size(T, 2));
        y = -(p(1)*x + p(3))/p(2);
        
        line = plot(x, y);
        scatter(t(i,1), t(i, 2), 30, get(line, 'Color'), 'filled')
        
        
    end
    
    hold off
    drawnow

end

