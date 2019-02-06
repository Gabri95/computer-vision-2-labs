function [r] = decompose_rotation(R)
  %{	
    x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
    r = [x, y, z];
%}
    sy = sqrt(R(1,1) * R(1,1) +  R(2,1) * R(2,1));
     
 
    if  sy >= 1e-6
        x = atan2(R(3,2) , R(3,3));
        y = atan2(-R(3,1), sy);
        z = atan2(R(2,1), R(1,1));
    else
        x = atan2(-R(2,3), R(2,2));
        y = atan2(-R(3,1), sy);
        z = 0;
    end
    r = [x, y, z];
end