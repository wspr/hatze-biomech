function R = rotation_matrix_zyx(r)
% ZYX rotation matrix: first around X, then Y, then Z

cx = cosd(r(1)); cy = cosd(r(2)); cz = cosd(r(3));
sx = sind(r(1)); sy = sind(r(2)); sz = sind(r(3));
R = ...
    [cz -sz 0;sz cz 0; 0 0 1]*...
    [cy 0 sy;0 1 0;-sy 0 cy]*...
    [1 0 0; 0 cx -sx; 0 sx cx];

end
