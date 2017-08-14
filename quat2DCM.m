function result = quat2DCM(q)
    if (size(q,1)==1 && size(q,2)==4)||(size(q,1)==4 && size(q,2)==1)
       a=q(1);b=q(2);c=q(3);d=q(4);
       result = [a^2+b^2-c^2-d^2, 2*(b*c-a*d), 2*(b*d+a*c);
                 2*(b*c+a*d), a^2-b^2+c^2-d^2, 2*(c*d-a*b);
                 2*(b*d-a*c), 2*(c*d+a*b), a^2-b^2-c^2+d^2;];
    else
       error('put in quaternion'); 
    end


end