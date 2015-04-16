function Ft = JacobianF( qnb, vn, pos, fb, xk )
%JACOB Summary of this function goes here
%   Detailed explanation goes here
    global glv
    [wnie, wnen, gn, retp] = earth(pos, vn);
    sl = retp.sl; cl = retp.cl; tl = retp.tl; 
    secl = 1./retp.cl; secl2 = secl^2;
    f_RMh = 1./retp.rmh; f_RNh = 1./retp.rnh;
    f_RMh2= f_RMh^2; f_RNh2 = f_RNh^2;
    Cnb = q2cnb(qnb);
    fn = qmulv(qnb,fb);
    M11 = [                          0,                                             glv.wie*sl+vn(1)*tl*f_RNh,                              -sin(xk(3))*vn(2)*f_RMh-cos(xk(3))*(glv.wie*cl+vn(1)*f_RNh);
                        -(glv.wie*sl+vn(1)*tl*f_RNh),                                           0,                                              -vn(2)*f_RMh+sin(xk(3))*(glv.wie*cl+vn(1)*f_RNh);
            sin(xk(3))*vn(2)*f_RMh+cos(xk(3))*(glv.wie*cl+vn(1)*f_RNh), cos(xk(3))*vn(2)*f_RMh-sin(xk(3))*(glv.wie*cl+vn(1)*f_RNh), (-xk(2)*sin(xk(3))+xk(1)*cos(xk(3)))*vn(2)*f_RMh-(xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)]; 
    M12 = [    0,    -f_RMh,    0;
             f_RNh,    0,       0; 
           tl*f_RNh,   0,       0];
    M13 = [           0,                 0,    vn(2)*f_RMh2;
                -glv.wie*sl,             0,   -vn(1)*f_RNh2; 
          glv.wie*cl+vn(1)*secl2*f_RNh,  0,  -vn(1)*tl*f_RNh2];
    M14 = -Cnb;
    M15 = zeros(3,3);
    M21 = [             0,                                      -fn(3),                                     cos(xk(3))*fn(2)-sin(xk(3))*fn(1);
                        fn(3),                                      0,                                      -sin(xk(3))*fn(2)-cos(xk(3))*fn(1);
            sin(xk(3))*fn(1)-cos(xk(3))*fn(2),      cos(xk(3))*fn(1)+sin(xk(3))*fn(2),  (-xk(2)*sin(xk(3))+xk(1)*cos(xk(3)))*fn(1)+(xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*fn(2) ];
    M22 = -askew(2*wnie+wnen)+askew(vn)*M12;
    M23 = askew(vn)*[           0,                   0,    vn(2)*f_RMh2;
                          -2*glv.wie*sl,             0,   -vn(1)*f_RNh2; 
                    2*glv.wie*cl+vn(1)*secl2*f_RNh,  0,  -vn(1)*tl*f_RNh2];
    M24 = zeros(3,3);
    M25 = Cnb;
    M31 = zeros(3,3);
    M32 = [     0,      f_RMh, 0;
            secl*f_RNh,   0,   0; 
                0,        0,   1];
    M33 = [           0,           0,       -vn(2)*f_RMh2; 
            vn(1)*tl*secl*f_RNh,   0,     -vn(1)*secl*f_RNh2; 
                      0,           0,             0];
    M34 = zeros(3,3);
    M35 = zeros(3,3);
    Melse = zeros(6,15);
    
    Ft = [ M11, M12, M13, M14, M15;
           M21, M22, M23, M24, M25;
           M31, M32, M33, M34, M35;
           Melse                   ];

end

