function result = EulerAccel(ax, ay, g)
    theta = asin(ax/g);
    phi = asin(-ay/(g*cos(theta)));

    result = [theta, phi];
end
