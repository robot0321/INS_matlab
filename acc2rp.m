function result = acc2rp(acc)
    if(size(acc,1)==3)
        ax = acc(1,:);
        ay = acc(2,:);
        az = acc(3,:);
    elseif(size(acc,2)==3)
        ax = acc(:,1);
        ay = acc(:,2);
        az = acc(:,3);
    else
        error('put in 3d accel vector');
    end
        
    theta = atan(ax/norm(ay,az));
    phi = atan2(-ay,-az);
    
    
    result = [phi, theta];
    
end
