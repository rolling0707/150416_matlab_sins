function [x, y, z] = DCMtoEuler(Cnb)%(秦永元P297--9.2.41)
%%%这是在方位角的正方向为北偏西的情况下进行推导的%%%
    x = asin(Cnb(3,2));
    
    y1 = atan(-Cnb(3,1)/Cnb(3,3));
    if Cnb(3,3)>0
        y = y1;
    else
        if y1>0
            y = y1-pi;
        else
            y = y1+pi;
        end
    end
    
    z1 = atan(-Cnb(1,2)/Cnb(2,2));
    if abs(Cnb(2,2))<1.e-20;
        if Cnb(1,2)>0
            z = -pi/2;
        else
            z = pi/2;
        end
    elseif Cnb(2,2)>0
        z = z1;
    else
        if Cnb(1,2)>0
            z = (z1-pi);
        else
            z = (z1+pi);
        end
    end
            
     