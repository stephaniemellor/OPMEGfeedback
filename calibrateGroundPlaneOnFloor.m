function [T,R] = calibrateGroundPlaneOnFloor(groundMarkers)
    
    %Calculate lengths of triangle lengths to work out which markers are which sides
    d12 = norm(groundMarkers(:,1) - groundMarkers(:,2));
    d23 = norm(groundMarkers(:,2) - groundMarkers(:,3));
    d31 = norm(groundMarkers(:,3) - groundMarkers(:,1));
    
    % Choose origin point and z and y sides
    if d12 < d31 && d23 < d31
        origin_idx = 2;
        if d12 < d23
            x_idx = 1;
            z_idx = 3;
        else
            x_idx = 3;
            z_idx = 1;
        end
    elseif d12 < d23 && d31 < d23
        origin_idx = 1;
        if d12 < d31
            x_idx = 2;
            z_idx = 3;
        else
            x_idx = 3;
            z_idx = 2;
        end
    else
        origin_idx = 3;
        if d23 < d31
            x_idx = 2;
            z_idx = 1;
        else
            x_idx = 1;
            z_idx = 2;
        end
    end
    
    % Move all points to the new origin point
    T = reshape(-groundMarkers(:,origin_idx),3,1);
    groundMarkers = groundMarkers + repmat(T,1,3);
    
    
    % Rotate frame so z lies along long side and x along shorter side
    % Rotate around z until z marker lies on x = 0
    sintheta = groundMarkers(1,z_idx)/sqrt(groundMarkers(1,z_idx)^2 + groundMarkers(2,z_idx)^2);
    costheta = groundMarkers(2,z_idx)/sqrt(groundMarkers(1,z_idx)^2 + groundMarkers(2,z_idx)^2);
    Rz = [costheta, -sintheta, 0; sintheta, costheta, 0; 0, 0, 1];
    groundMarkers = Rz*groundMarkers;
    
    % Rotate around x so z marker lies on y = 0
    sintheta = groundMarkers(2,z_idx)/sqrt(groundMarkers(3,z_idx)^2 + groundMarkers(2,z_idx)^2);
    costheta = groundMarkers(3,z_idx)/sqrt(groundMarkers(3,z_idx)^2 + groundMarkers(2,z_idx)^2);
    Rx = [1, 0, 0; 0, costheta, -sintheta; 0, sintheta, costheta];
    groundMarkers = Rx*groundMarkers;

    % Create a y marker to ensure coordinate system is right-handed
    ymarker = cross(groundMarkers(:,z_idx), groundMarkers(:,x_idx));
    
    % Rotate around z until y marker lies on x = 0
    sintheta = ymarker(1)/sqrt(ymarker(1)^2 + ymarker(2)^2);
    costheta = ymarker(2)/sqrt(ymarker(1)^2 + ymarker(2)^2);
    Rz1 = [costheta, -sintheta, 0; sintheta, costheta, 0; 0, 0, 1];
    
    % Create overall matrix
    R = Rz1*Rx*Rz;

    % Update Translation
    T = R*T;
    
    end