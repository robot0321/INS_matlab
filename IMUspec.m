function [a,b,c,d]=IMUspec(name)
    if name=='HG1700AG58'
        a = 1;
        b = 0.125;
        c = 1;
        d = 0.65;
    else
        error('wrong name');
    end
end