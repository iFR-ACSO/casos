function [q] = EulerAngle2Quat321(roll,pitch,yaw)

        q = [cos(roll/2)*cos(pitch/2)*cos(yaw/2) + sin(roll/2)*sin(pitch/2)*sin(yaw/2);
            -cos(roll/2)*sin(pitch/2)*sin(yaw/2) + cos(pitch/2)*cos(yaw/2)*sin(roll/2);
             cos(roll/2)*cos(yaw/2)*sin(pitch/2) + sin(roll/2)*cos(pitch/2)*sin(yaw/2);
             cos(roll/2)*cos(pitch/2)*sin(yaw/2) - sin(roll/2)*cos(yaw/2)*sin(pitch/2)];

end

