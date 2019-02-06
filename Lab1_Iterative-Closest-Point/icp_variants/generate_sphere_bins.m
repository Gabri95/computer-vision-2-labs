

function vectors = generate_sphere_bins(N)
%     https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
%     vectors = zeros(3, N);
% 
%     count = 0;
%     a = 4 * pi/N;
%     d = sqrt(a);
% 
%     M_theta = round(pi/d);
%     d_theta = pi/M_theta;
%     d_phi = a/ d_theta;
%     for m = 0:M_theta-1
%         theta = pi*(m + 0.5)/M_theta;
%         M_phi = round(2*pi*sin(theta/d_phi));
%         for n =0:M_phi - 1
%             phi = 2*pi*n / M_phi;
%             vectors(:, count + 1) = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)];
%             count = count + 1;
%         end
%     end
%     count


    vectors = zeros(3, N);
    offset = 2./N;
    increment = pi * (3 - sqrt(5));
    rnd = 1; %randi(N);
    
    for i = 0:N-1
        y = ((i * offset) - 1) + (offset / 2);
        r = sqrt(1 - y^2);

        phi = mod(i + rnd, N) * increment;

        x = cos(phi) * r;
        z = sin(phi) * r;

        vectors(:, i+1) = [x;y;z];
    end
end

