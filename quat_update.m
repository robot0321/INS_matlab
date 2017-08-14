function new_q = quat_update(q,w,dt)
    if (size(q,1)==1 && size(q,2)==4)||(size(q,1)==4 && size(q,2)==1)
        if (size(w,1)==1 && size(w,2)==3)||(size(w,1)==3 && size(w,2)==1)
        a=q(1);b=q(2);c=q(3);d=q(4);
        wx=w(1);wy=w(2);wz=w(3);
        % sigma completion
        
        % q
%         sigma = w*dt;
%         
%         ac = cos(norm(sigma)/2);
%         if norm(sigma)==0
%             as = 0;
%         else
%             as = sin(norm(sigma)/2)/norm(sigma);
%         end
%         r = [ac, as*sigma(1), as*sigma(2), as*sigma(3)]';
%         new_q = product_q(q,r);

        delta_q=zeros(4,1);
        delta_q(1) = -0.5*(b*wx + c*wy + d*wz);
        delta_q(2) = 0.5*(a*wx - d*wy + c*wz);
        delta_q(3) = 0.5*(d*wx + a*wy - b*wz);
        delta_q(4) = -0.5*(c*wx - b*wy - a*wz);

        new_q = q + delta_q*dt;
        else
            error('put in angular velocity'); 
        end
    else
       error('put in quaternion'); 
    end
end