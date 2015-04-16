tic;
clear all
clc
glvs
Ts = 0.005;
imu0 = load('data\RAWIMU.txt');
veli = load('data\BESTVEL.txt');
posi = load('data\BESTPOS.txt');
pos0 = [40.067987*glv.deg; 116.247871*glv.deg; 43.53];
vn0 = [0; 0; 0];
[pitch, roll, yaw] = coarse_alignment(imu0(1:47000, :), pos0);
att0 = [pitch; roll; yaw];
cnb0 = a2cnb(att0);
qnb0 = m2qnb(cnb0);
qnb = qnb0; vn = vn0; pos = pos0;
n = 1;                                                  %n子样数
Ts = n*Ts; 

dpos = [1./glv.Re; 1./glv.Re; .3];
dv = [.5; .5; .1];
[eb, web, db, wdb] = drift_bias([1;1;1;1]);
%%%%滤波器参数设置%%%%
Qt = diag([web; wdb; repmat(0.01, 3, 1); 0.1*web; 0.1*wdb])^2;
Rk = diag([dv; dpos])^2;
Pk = diag([[10.;10.;30.]*glv.deg; [.2; .2; .1]; [5./glv.Re;...
    5./glv.Re; 2.]; [2.; 2.; 2.]*glv.dph; ...
    [3.; 3.; 3.]*glv.mg])^2;
Xk = [[5.;5.;25.]*glv.deg; [.3; .3; .1]; ...
    [5./glv.Re; 5./glv.Re; 1.]; [.7; .7; .7]*glv.dph;...
    [1.; 1.; 1.]*glv.mg];
Hk = [zeros(6, 3), eye(6), zeros(6,6)];
len = 252760;
lenp = 1020;
Xkk_ukf = zeros(lenp, length(Xk));
Zk_ukf = zeros(lenp, 6);
navigation = zeros(len, 9);
kv = 1;
kk = 1;
kerr = 1;
for k=47000:n:299750
    k1 = k+n-1;
    wm = imu0(k:k1, 6:8)';
    vm = imu0(k:k1, 3:5)';
   [qnb, vn, pos] = sinss(qnb, vn, pos, wm, vm, Ts);
   cnb = q2cnb(qnb);
   [x, y, z] = DCMtoEuler(cnb);
   att = [x; y; z];
   navigation(kk,:) = [att'./glv.deg, vn', pos(1:2,1)'/glv.deg, pos(3,1)];
   wm = sum(wm,2);
   vm = sum(vm,2);
   fb = vm./Ts;
   [wnie, wnen, gn, retp] = earth(pos, vn);
   sl = retp.sl; cl = retp.cl; tl = retp.tl;
   secl = 1./retp.cl; secl2 = secl^2;
   f_RMh = 1./retp.rmh; f_RNh = 1./retp.rnh;
   f_RMh2 = f_RMh^2; f_RNh2 = f_RNh^2;
   fn = qmulv(qnb, fb);
   Cnb = q2cnb(qnb);
   T = Ts;
%    fstate = @(xk)[xk(1)+((cos(xk(3))-1)*vn(2)*f_RMh-sin(xk(3))*(glv.wie*cl+vn(1)*f_RNh)+...
%         xk(2)*(glv.wie*sl+vn(1)*tl*f_RNh)-f_RMh*xk(5)+vn(2)*f_RMh2*xk(9)+...
%         -cos(xk(3))*(xk(10)*Cnb(1,1)+xk(11)*Cnb(1,2)+xk(12)*Cnb(1,3))-...
%         sin(xk(3))*(xk(10)*Cnb(2,1)+xk(11)*Cnb(2,2)+xk(12)*Cnb(2,3)))*T;
%         xk(2)+(-xk(3)*vn(2)*f_RMh+(1-cos(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)+...
%         -xk(1)*(glv.wie*sl+vn(1)*tl*f_RNh)-glv.wie*sl*xk(7)+f_RNh*xk(4)+...
%         -vn(1)*f_RNh2*xk(9)+sin(xk(3))*(xk(10)*Cnb(2,1)+xk(11)*Cnb(2,2)+xk(12)*Cnb(2,3))...
%         -cos(xk(3))*xk(10)*Cnb(2,1)+xk(11)*Cnb(2,2)+xk(12)*Cnb(2,3))*T;
%         xk(3)+((xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*vn(2)*f_RMh+...
%         (-xk(2)*sin(xk(3))+xk(1)*cos(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)+...
%         (glv.wie*cl+vn(1)*secl2*f_RNh)*xk(7)+tl*f_RNh*xk(4)+...
%         -vn(1)*tl*f_RNh2*xk(9)-(xk(10)*Cnb(3,1)+xk(11)*Cnb(3,2)+xk(12)*Cnb(3,3)))*T;
%         xk(4)+((cos(xk(3))-1)*fn(1)+sin(xk(3))*fn(2)-xk(2)*fn(3)+(vn(2)*tl-vn(3))*f_RNh*xk(4)+...
%         (2*glv.wie*sl+vn(1)*tl*f_RNh)*xk(5)-(2*glv.wie*cl+vn(1)*f_RNh)*xk(6)+...
%         (2*glv.wie*(vn(3)*sl+vn(2)*cl)+vn(1)*vn(2)*secl2*f_RNh)*xk(7)+...
%         (vn(1)*vn(3)-vn(1)*vn(2)*tl)*f_RNh2*xk(9)+cos(xk(3))*(xk(13)*Cnb(1,1)+xk(14)*Cnb(1,2)+xk(15)*Cnb(1,3))+...
%         sin(xk(3))*(xk(13)*Cnb(2,1)+xk(14)*Cnb(2,2)+xk(15)*Cnb(2,3)))*T;
%         xk(5)+(-sin(xk(3))*fn(1)+(cos(xk(1))-1)*fn(2)+xk(1)*fn(3)-(2*glv.wie*sl+vn(1)*tl*f_RNh)*xk(4)+...
%         -vn(3)*f_RMh*xk(5)-vn(2)*f_RMh*xk(6)-(2*glv.wie*vn(1)*cl+vn(1)^2*secl2*f_RNh)*xk(7)+...
%         (vn(2)*vn(3)*f_RMh2+vn(1)^2*tl*f_RNh2)*xk(9)+cos(xk(3))*(xk(13)*Cnb(2,1)+xk(14)*Cnb(2,2)+xk(15)*Cnb(2,3))...
%         -sin(xk(3))*(xk(13)*Cnb(1,1)+xk(14)*Cnb(1,2)+xk(15)*Cnb(1,3)))*T;
%         xk(6)+((xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*fn(1)+(xk(2)*sin(xk(3))-xk(1)*cos(xk(3)))*fn(2)+...
%         2*(glv.wie*cl+vn(1)*f_RNh)*xk(4)+2*vn(2)*f_RMh*xk(5)-2*glv.wie*vn(1)*sl*xk(7)+...
%         -(vn(2)^2*f_RMh2+vn(1)^2*f_RNh2)*xk(9)+(xk(13)*Cnb(3,1)+xk(14)*Cnb(3,2)+xk(15)*Cnb(3,3)))*T;
%         xk(7)+(xk(5)*f_RMh-xk(9)*vn(2)*f_RMh2)*T;
%         xk(8)+(xk(4)*secl*f_RNh+xk(7)*vn(1)*tl*secl*f_RNh-xk(9)*vn(1)*secl*f_RNh2)*T;
%         xk(9)+(xk(6))*T;
%         xk(10);xk(11);xk(12);xk(13);xk(14);xk(15)];
   fstate1 = @(xk)[((cos(xk(3))-1)*vn(2)*f_RMh-sin(xk(3))*(glv.wie*cl+vn(1)*f_RNh)+...
        xk(2)*(glv.wie*sl+vn(1)*tl*f_RNh)-f_RMh*xk(5)+vn(2)*f_RMh2*xk(9)+...
        -cos(xk(3))*(xk(10)*Cnb(1,1)+xk(11)*Cnb(1,2)+xk(12)*Cnb(1,3))-...
        sin(xk(3))*(xk(10)*Cnb(2,1)+xk(11)*Cnb(2,2)+xk(12)*Cnb(2,3)));
        (-xk(3)*vn(2)*f_RMh+(1-cos(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)+...
        -xk(1)*(glv.wie*sl+vn(1)*tl*f_RNh)-glv.wie*sl*xk(7)+f_RNh*xk(4)+...
        -vn(1)*f_RNh2*xk(9)+sin(xk(3))*(xk(10)*Cnb(2,1)+xk(11)*Cnb(2,2)+xk(12)*Cnb(2,3))...
        -cos(xk(3))*xk(10)*Cnb(2,1)+xk(11)*Cnb(2,2)+xk(12)*Cnb(2,3));
        ((xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*vn(2)*f_RMh+...
        (-xk(2)*sin(xk(3))+xk(1)*cos(xk(3)))*(glv.wie*cl+vn(1)*f_RNh)+...
        (glv.wie*cl+vn(1)*secl2*f_RNh)*xk(7)+tl*f_RNh*xk(4)+...
        -vn(1)*tl*f_RNh2*xk(9)-(xk(10)*Cnb(3,1)+xk(11)*Cnb(3,2)+xk(12)*Cnb(3,3)));
        ((cos(xk(3))-1)*fn(1)+sin(xk(3))*fn(2)-xk(2)*fn(3)+(vn(2)*tl-vn(3))*f_RNh*xk(4)+...
        (2*glv.wie*sl+vn(1)*tl*f_RNh)*xk(5)-(2*glv.wie*cl+vn(1)*f_RNh)*xk(6)+...
        (2*glv.wie*(vn(3)*sl+vn(2)*cl)+vn(1)*vn(2)*secl2*f_RNh)*xk(7)+...
        (vn(1)*vn(3)-vn(1)*vn(2)*tl)*f_RNh2*xk(9)+cos(xk(3))*(xk(13)*Cnb(1,1)+xk(14)*Cnb(1,2)+xk(15)*Cnb(1,3))+...
        sin(xk(3))*(xk(13)*Cnb(2,1)+xk(14)*Cnb(2,2)+xk(15)*Cnb(2,3)));
        (-sin(xk(3))*fn(1)+(cos(xk(1))-1)*fn(2)+xk(1)*fn(3)-(2*glv.wie*sl+vn(1)*tl*f_RNh)*xk(4)+...
        -vn(3)*f_RMh*xk(5)-vn(2)*f_RMh*xk(6)-(2*glv.wie*vn(1)*cl+vn(1)^2*secl2*f_RNh)*xk(7)+...
        (vn(2)*vn(3)*f_RMh2+vn(1)^2*tl*f_RNh2)*xk(9)+cos(xk(3))*(xk(13)*Cnb(2,1)+xk(14)*Cnb(2,2)+xk(15)*Cnb(2,3))...
        -sin(xk(3))*(xk(13)*Cnb(1,1)+xk(14)*Cnb(1,2)+xk(15)*Cnb(1,3)));
        ((xk(2)*cos(xk(3))+xk(1)*sin(xk(3)))*fn(1)+(xk(2)*sin(xk(3))-xk(1)*cos(xk(3)))*fn(2)+...
        2*(glv.wie*cl+vn(1)*f_RNh)*xk(4)+2*vn(2)*f_RMh*xk(5)-2*glv.wie*vn(1)*sl*xk(7)+...
        -(vn(2)^2*f_RMh2+vn(1)^2*f_RNh2)*xk(9)+(xk(13)*Cnb(3,1)+xk(14)*Cnb(3,2)+xk(15)*Cnb(3,3)));
        (xk(5)*f_RMh-xk(9)*vn(2)*f_RMh2);
        (xk(4)*secl*f_RNh+xk(7)*vn(1)*tl*secl*f_RNh-xk(9)*vn(1)*secl*f_RNh2);
        (xk(6));
        0;0;0;0;0;0];
    fstate = discretize(fstate1, Ts);
    hmeas = @(xk)[xk(4);xk(5);xk(6);xk(7);xk(8);xk(9)];
    Qk = Qt*Ts;
   for j = k:1:k1
       if imu0(j, 2)==veli(kv, 2)
           vi = veli(kv, 3:5)';
           pi = [posi(kv, 3:4)'*glv.deg; posi(kv, 5)];
           Zk = [vn-vi; pos-pi];
           Zk_ukf(kerr, :) = Zk';
           [Xk, Pk] = ukf(fstate, Qk, Xk, Pk, hmeas, Zk, Rk);
           Xkk_ukf(kerr, :) = Xk';
           vn = vn-Xk(4:6, 1);
           pos = pos-Xk(7:9, 1);
           pos(3) = posi(kv, 5);
           xk(9)=0;
           kv = kv+1;
           kerr = kerr+1;
       else
           [Xk, Pk] = ukf(fstate, Qk, Xk, Pk);
       end
   end
   timedisp(kk, Ts, 100)
   kk = kk+1;
end

save alignment_ukf_rk4 navigation
save Xkk_ukf3_rk4 Xkk_ukf
save Zk_ukf_rk4 Zk_ukf

figure,             %画出组合导航轨迹
plot(navigation(1:kk-1,8), navigation(1:kk-1,7), '-'), grid on
xlabel('\it\lambda\rm / \circ'); ylabel('\itLatitude\rm / \circ');

figure,             %画出俯仰误差角曲线
plot(Xkk_ukf(1:kerr-1,1)/glv.deg, '-b'), grid on
xlabel('\itt/s'); ylabel('\it\theta\rm / \circ'); grid on

figure,             %画出俯横滚差角曲线
plot(Xkk_ukf(1:kerr-1,2)/glv.deg, '-b'), grid on
xlabel('\itt/s'); ylabel('\it\gamma\rm / \circ');

figure,             %画出俯方位差角曲线
plot(Xkk_ukf(1:kerr-1,3)/glv.deg, '-b'), grid on
xlabel('\itt/s'); ylabel('\it\psi\rm / \circ');

figure,             %画出东北向速度误差
subplot(2,1,1), plot(Zk_ukf(1:kerr-1,1)), grid on
xlabel('\itt/s'); ylabel('\it\delta Ve\rm /m/s'); 
subplot(2,1,2), plot(Zk_ukf(1:kerr-1,2)), grid on
xlabel('\itt/s'); ylabel('\it\delta Vn\rm /m/s'); 

figure,             %画出纬度、经度误差
subplot(2,1,1), plot(Zk_ukf(1:kerr-1,4)*glv.Re), grid on
xlabel('\itt/s'); ylabel('\it\delta Latitude\rm /m'); 
subplot(2,1,2), plot(Zk_ukf(1:kerr-1,5)*glv.Re), grid on
xlabel('\itt/s'); ylabel('\it\delta lambda\rm /m'); 

toc
