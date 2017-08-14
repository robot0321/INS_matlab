function result = acc2att(acc)
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
        
    abmag = sqrt(ax.^2 + ay.^2 + az.^2);
    ax = ax./abmag;
    ay = ay./abmag;
    az = az./abmag;

    % phi = atan(ay./az);
    % theta = asin(-ax/1); % g ¥‹¿ß
    % ly = mx.*cos(theta)*sin(p%my.*cos(phi) - mz.*sin(phi);
    % lz = -mx.*sin(theta) + my.*cos(theta).*sin(phi) + mz.*cos(theta).*cos(phi);
    % psi = atan(ly./lz);
    phi = atan(ay./sqrt(ax.^2+az.^2));
    theta = atan(ax./sqrt(ay.^2+az.^2));
    psi = atan(az./sqrt(ax.^2+az.^2));


    for i=2:size(phi)
       if(abs(phi(i)+pi - phi(i-1)) < abs(phi(i)-phi(i-1)))
          phi(i) = phi(i) + pi;
       end
       if(abs(phi(i)-pi - phi(i-1)) < abs(phi(i)-phi(i-1)))
           phi(i) = phi(i) - pi;
       end
    end

    for i=2:size(psi)
       if(abs(psi(i)+pi - psi(i-1)) < abs(psi(i)-psi(i-1)))
          psi(i) = psi(i) + pi;
       end
       if(abs(psi(i)-pi - psi(i-1)) < abs(psi(i)-psi(i-1)))
           psi(i) = psi(i) - pi;
       end
    end
    
    result = [phi, theta, psi];
end
