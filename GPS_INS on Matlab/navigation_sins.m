tic;
clear all
clc
glvs
Ts = 0.005;
imu0 = load('data\RAWIMU.txt');
veli = load('data\BESTVEL.txt');
posi = load('data\BESTPOS.txt');
pos0 = [40.0678531*glv.deg; 116.247326*glv.deg; 47.92];
vn0 = [0; 0; 0];
att0 = [0.0081; -0.0065; -2.9538];
% att0 = [0.0040; -0.0048; -2.9976]; %采用UKF滤波对准最后阶段的姿态值
cnb0 = a2cnb(att0);
qnb0 = m2qnb(cnb0);
qnb = qnb0; vn = vn0; pos = pos0;
n = 4;
Ts = n*Ts;

len = length(imu0);
lenp = length(veli);
Zk_sins = zeros(lenp, 6);
navigation = zeros(len, 9);
kv = 2035;
kk = 1;
kerr = 1;
for k = 454321:n:716560
   k1 = k+n-1;
   wm = imu0(k:k1,6:8)';
   vm = imu0(k:k1,3:5)';
   [qnb, vn, pos] = sinss(qnb, vn, pos, wm, vm, Ts);
   cnb = q2cnb(qnb);
   [x, y, z] = DCMtoEuler(cnb);
   att = [x; y; z];
   navigation(kk,:) = [att'./glv.deg, vn', pos(1:2,1)'/glv.deg, pos(3,1)];
   for j=k:1:k1
           if imu0(j, 2)==veli(kv, 2)               
           vi = veli(kv, 3:5)';
           pi = [posi(kv, 3:4)'*glv.deg; posi(kv, 5)];
           Zk = [vn-vi; pos-pi];
           Zk_sins(kerr, :) = Zk';
           kv = kv+1;
           kerr = kerr+1;
           end
   end
   kk = kk+1;
   timedisp(k, Ts, 100);
end

save navigation_sins navigation

save Zk_sins Zk_sins

figure,             %画出惯性导航轨迹
plot(navigation(1:kk-1,8), navigation(1:kk-1,7), '-'),
xlabel('\it\lambda\rm / \circ'); ylabel('\itLatitude\rm / \circ');

figure,             %画出东北向速度误差
subplot(2,1,1), plot(Zk_sins(1:kerr-1,1),'-'),
xlabel('\itt/s'); ylabel('\it\delta Ve\rm /m/s'); grid on
subplot(2,1,2), plot(Zk_sins(1:kerr-1,2), '-'),
xlabel('\itt/s'); ylabel('\it\delta Vn\rm /m/s'); grid on

figure,             %画出纬度、经度误差
subplot(2,1,1), plot(Zk_sins(1:kerr-1,4)*glv.Re, '-'),
xlabel('\itt/s'); ylabel('\it\delta Latitude\rm /m'); grid on
subplot(2,1,2), plot(Zk_sins(1:kerr-1,5)*glv.Re, '-'),
xlabel('\itt/s'); ylabel('\it\delta lambda\rm /m'); grid on

toc


