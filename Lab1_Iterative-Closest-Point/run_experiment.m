function error = run_experiment(folder,source, target, source_normals, target_normals, sampling, selection, matching, rejection, weighting, plot, real_rot, real_trans)

    fprintf('Sampling: %s; Selection: %s; Matching: %s; Rejection: %s; Weighting: %s\n', sampling{1}, selection{1}, matching{1}, rejection{1}, weighting{1});
    tic;
    [transformed, R, t, error, iterations] = icp(source, target, source_normals, target_normals, sampling{2}, selection{2}, matching{2}, rejection{2}, weighting{2}, plot);
    elapsedTime = toc;
    if nargin > 10 
        translation_error = sum(abs(real_trans - t'));
        rotation_error = sum(abs(decompose_rotation(R) - real_rot));
        fprintf("Rotation: [%d %d %d]: [%d %d %d] \n" ,decompose_rotation(R),real_rot);
        fprintf("Translation: [%d %d %d]: [%d %d %d] \n" ,t',real_trans);
        
        fprintf("Rotation error: %d \n" ,rotation_error);
        fprintf("Translation error: %d \n" ,translation_error);
    end
    
    fig = figure;
    
    plotPCL(R*source + t, [])
    hold on
    plotPCL(target, [])
    
    plot_title = sprintf('Final Error: %f', error);
    plot_subtitle = sprintf('Total Iterations: %d  Time: %d', iterations,elapsedTime);
    
    title({plot_title, plot_subtitle})
    
    hold off 
    
    
    name = fullfile(folder,sprintf('%s_%s_%s_%s_%s', sampling{1}, selection{1}, matching{1}, rejection{1}, weighting{1}));
    saveas(fig,sprintf('%s.png',name));
    close all
    
    if nargin > 10        
        filename = sprintf('%s.txt', name);
        
        fileID = fopen(filename, 'w');
        fprintf(fileID, 'Translation Error: %0.5f\n', translation_error);
        fprintf(fileID, 'Rotation Error: %0.5f\n', rotation_error);
        %fprintf(fileID, 'Rotation Angles: %0.5f, %0.5f, %0.5f\n', rotation_error_angles(1), rotation_error_angles(2), rotation_error_angles(3));
        
        fclose(fileID);
        
    end
    
end
