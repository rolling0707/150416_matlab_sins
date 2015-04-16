function fx = efff( qnb, vn, pos, fb, xk )
%%%得到导数值
    global glv
    [wnie, wnen, gn, retp] = earth(pos, vn);
    sl = retp.sl; cl = retp.cl; tl = retp.tl; 
    secl = 1./retp.cl; secl2 = secl^2;
    f_RMh = 1./retp.rmh; f_RNh = 1./retp.rnh;
    f_RMh2= f_RMh^2; f_RNh2 = f_RNh^2;
    Cnb = q2cnb(qnb);
    xk(10:12) = qmulv(qnb, xk(10:12));
    xk(13:15) = qmulv(qnb, xk(13:15));
    fn = qmulv(qnb,fb);
    m1 = (cos(xk(3))-1)*vn(2)*f_RMh-sin(xk(3))*(glv.wie*cl+vn(1)*f_RNh)+...
        xk(2)*(glv.wie*sl+vn(1)*tl*f_RNh)-f_RMh*xk(5)+vn(2)*f_RMh2*xk(9)+...
        -cos(xk(3))*xk(10)-sin(xk(3))*xk(11);
    %%
    m2 = -xk(3)*vn(2)*f_RMh+(1-cos(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)+...
        -xk(1)*(glv.wie*sl+vn(1)*tl*f_RNh)-glv.wie*sl*xk(7)+f_RNh*xk(4)+...
        -vn(1)*f_RNh2*xk(9)+sin(xk(3))*xk(10)-cos(xk(3))*xk(11);
    %%
    m3 = (xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*vn(2)*f_RMh+...
        (-xk(2)*sin(xk(3))+xk(1)*cos(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)+...
        (glv.wie*cl+vn(1)*secl2*f_RNh)*xk(7)+tl*f_RNh*xk(4)+...
        -vn(1)*tl*f_RNh2*xk(9)-xk(12);
    %%
    m4 = (cos(xk(3))-1)*fn(1)+sin(xk(3))*fn(2)-xk(2)*fn(3)+(vn(2)*tl-vn(3))*f_RNh*xk(4)+...
        (2*glv.wie*sl+vn(1)*tl*f_RNh)*xk(5)-(2*glv.wie*cl+vn(1)*f_RNh)*xk(6)+...
        (2*glv.wie*(vn(3)*sl+vn(2)*cl)+vn(1)*vn(2)*secl2*f_RNh)*xk(7)+...
        (vn(1)*vn(3)-vn(1)*vn(2)*tl)*f_RNh2*xk(9)+cos(xk(3))*xk(13)+sin(xk(3))*xk(14);
    %%
    m5 = -sin(xk(3))*fn(1)+(cos(xk(1))-1)*fn(2)+xk(1)*fn(3)-(2*glv.wie*sl+vn(1)*tl*f_RNh)*xk(4)+...
        -vn(3)*f_RMh*xk(5)-vn(2)*f_RMh*xk(6)-(2*glv.wie*vn(1)*cl+vn(1)^2*secl2*f_RNh)*xk(7)+...
        (vn(2)*vn(3)*f_RMh2+vn(1)^2*tl*f_RNh2)*xk(9)+cos(xk(3))*xk(14)-sin(xk(3))*xk(13);
    %%
    m6 = (xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*fn(1)+(xk(2)*sin(xk(3))-xk(1)*cos(xk(3)))*fn(2)+...
        2*(glv.wie*cl+vn(1)*f_RNh)*xk(4)+2*vn(2)*f_RMh*xk(5)-2*glv.wie*vn(1)*sl*xk(7)+...
        -(vn(2)^2*f_RMh2+vn(1)^2*f_RNh2)*xk(9)+xk(15);
    %%
    m7 = xk(5)*f_RMh-xk(9)*vn(2)*f_RMh2;
    %%
    m8 = xk(4)*secl*f_RNh+xk(7)*vn(1)*tl*secl*f_RNh-xk(9)*vn(1)*secl*f_RNh2;
    %%
    m9 = xk(6);
    m10 = 0;
    m11 = 0;
    m12 = 0;
    m13 = 0;
    m14 = 0;
    m15 = 0;
    fx = [m1, m2, m3, m4, m5, m6, m7, m8, m9,...
        m10, m11, m12, m13, m14, m15]';


end

