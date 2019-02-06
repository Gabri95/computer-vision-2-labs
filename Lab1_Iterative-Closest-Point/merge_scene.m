
function [cloud, cloud_normals, error] = merge_scene(frames, frames_normals, rate, iterative_merge, reference, sampling, selection, matching, rejection, weighting, plot, name)
    
    F = length(frames);
    
    if isempty(reference)
        reference = 1 + rate* ceil(F/(2*rate));
    end
    
    f = 1;
    
    base = frames{f};
    base_normals = frames_normals{f};
    
    
    
    cloud = [base];
    cloud_normals = [base_normals];
    
    f = f + rate;
    
    
    FirstRotation = eye(3);
    FirstTranslation = zeros(3, 1);
    
    Rotation = eye(3);
    Translation = zeros(3, 1);
    Transformation = eye(4);
    
    
    if plot
        figure();
    end
    
    while f < F
        start_frame = tic;
        
        fprintf('Frame %d, base size: %d\n', f, size(base, 2));
        
        target = frames{f};
        target_normals = frames_normals{f};
        
        
        [transformed, R, t, ~, ~] = icp(base, target, base_normals, target_normals, sampling, selection, matching, rejection, weighting, false);
         
%         transf = eye(4);
%         transf(1:3, 1:3) = R;
%         transf(1:3, 4) = t;
%         Transformation = transf * Transformation; 
%         Inverse_Transformation = inv(Transformation);
%         R_1 = Inverse_Transformation(1:3, 1:3);
%         t_1 = Inverse_Transformation(1:3, 4);
%         transformed = R_1*target + t_1;
        
        
        if f <= reference
            
            FirstRotation = R*FirstRotation;
            FirstTranslation = R*FirstTranslation + t;
            
            cloud = [R*cloud + t, target];
            cloud_normals = [R*cloud_normals, target_normals];
            
        elseif f > reference
            
            Rotation = R*Rotation;
            Translation = R*Translation + t;
        
            transformed = inv(Rotation)*(target - Translation);
            cloud = [cloud, transformed];
            cloud_normals = [cloud_normals, inv(Rotation)*target_normals];
            
            
        end
        
        if plot
    %         plotPCL(cloud(:, (rand(1, size(cloud, 2)) - rate/(f-1)) <= 0), []);
    %         plotPCL(normalspace_samples(cloud, cloud_normals, rate/(f-1)), []);
    %         plotPCL(bucketnormalspace_samples(transformed, inv(Rotation)*target_normals, 0.01), []);
    
            transformed_base = R*base + t;        
            idxs_b = (rand(1, size(transformed_base, 2)) - 0.1) <= 0;
            idxs_t = (rand(1, size(target, 2)) - 0.1) <= 0;
            
            plotPCLs(transformed_base(:, idxs_b), target(:, idxs_t), false);
            drawnow
        end
        
        
        if iterative_merge
%             R = Transformation(1:3, 1:3);
%             t = Transformation(1:3, 4);
            
            p = 1 / ((f -1)/rate + 1);

            idxs = (rand(1, size(cloud, 2)) - p) <= 0;
            
            base = cloud(:, idxs);
            base_normals = cloud_normals(:, idxs);
            
            if f > reference
                base = Rotation*base + Translation;
                base_normals = Rotation*base_normals;
            end
            
%             plotPCL(base, []);
%             drawnow
        else
            base = target;
            base_normals = target_normals;
        end
        
        
        fprintf('Frame %d, time elapsed: %f\n', f, toc(start_frame));
        f = f + rate;
    end



    target = frames{1};
    target_normals = frames_normals{1};
    
    
    
    
    [transformed, R, t] = icp(base, target, base_normals, target_normals, sampling, selection, matching, rejection, weighting, false);
    
    figure()
    transformed_base = R*base + t;        
    idxs_b = (rand(1, size(transformed_base, 2)) - 0.1) <= 0;
    idxs_t = (rand(1, size(target, 2)) - 0.1) <= 0;
    plotPCLs(transformed_base(:, idxs_b), target(:, idxs_t), false);
    
    
    Rotation = R*Rotation;
    Translation = R*Translation + t;
    
    Rotation_identity = Rotation*FirstRotation
    Translation_identity = Rotation*FirstTranslation + Translation
    
    translation_error = sum(Translation_identity.^2)
    rotation_error = sum(sum((Rotation_identity - eye(3)).^2))
    
    error = rotation_error + translation_error;
    
    rotation_error_angles = decompose_rotation(Rotation_identity)
    
    back_transformed = inv(Rotation)*(target - Translation);
    forward_transformed = FirstRotation*target + FirstTranslation;
    
    if ~isempty(name) || plot
    
        fig=figure();
        plotPCL(cloud(:, (rand(1, size(cloud, 2)) - 0.01) <= 0), 'blue');
        hold on
        plotPCL(back_transformed(:, (rand(1, size(back_transformed, 2)) - 0.1) <= 0), 'red');
        plotPCL(forward_transformed(:, (rand(1, size(forward_transformed, 2)) - 0.1) <= 0), 'green');
        hold off
        
    end
    
    if ~isempty(name)
        
        saveas(fig, sprintf('%s.png', name));
        
        
        filename = sprintf('%s.txt', name);
        
        dlmwrite(filename, Rotation_identity);
        
        dlmwrite(filename, Translation_identity, '-append');
        
        fileID = fopen(filename, 'a');
        fprintf(fileID, 'Translation Error: %0.5f\n', translation_error);
        fprintf(fileID, 'Rotation Error: %0.5f\n', rotation_error);
        fprintf(fileID, 'Rotation Angles: %0.5f, %0.5f, %0.5f\n', rotation_error_angles(1), rotation_error_angles(2), rotation_error_angles(3));
        
        fclose(fileID);
        
    end

    if ~ plot
        close all
    end

end


function plotPCLs(A,B, new_figure)
    if new_figure
        figure;
    end
    plotPCL(A, 'red')
    hold on
    plotPCL(B, 'blue')
    hold off 
end

