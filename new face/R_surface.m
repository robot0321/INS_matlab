function result = R_surface(R_long, R_short, Latitude_rad)
    result = R_long*R_short/sqrt((R_short*cos(Latitude_rad)).^2 + (R_long*sin(Latitude_rad)).^2);
end
