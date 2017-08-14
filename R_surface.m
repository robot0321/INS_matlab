function result = R_surface(R_long, R_short, Latitude)
    result = R_long*R_short./sqrt((R_short*cos(Latitude)).^2 + (R_long*sin(Latitude)).^2);
end
