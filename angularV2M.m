function result = angularV2M(angular_Vector)
    if size(angular_Vector)==[3,1]
        result = [0, -angular_Vector(3), angular_Vector(2); angular_Vector(3), 0, -angular_Vector(1); -angular_Vector(2), angular_Vector(1), 0];
    elseif size(angular_Vector)==[1,3]
        result = [0, -angular_Vector(3), angular_Vector(2); angular_Vector(3), 0, -angular_Vector(1); -angular_Vector(2), angular_Vector(1), 0];
    else
        error('put in 3x1 or 1x3 angular vector');
    end   
end