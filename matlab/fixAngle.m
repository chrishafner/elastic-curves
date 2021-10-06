function ret = fixAngle(angle, angle_ref)
    ret = angle;
    while ret - angle_ref > pi
        ret = ret - 2*pi;
    end
    while ret - angle_ref < -pi
        ret = ret + 2*pi;
    end
end