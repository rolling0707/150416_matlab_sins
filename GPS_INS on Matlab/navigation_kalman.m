tic;
clear all
clc
glvs
Ts = 0.005;                                                 %IMU�����������
imu0 = load('data\RAWIMU.txt');
veli = load('data\BESTVEL.txt');
posi = load('data\BESTPOS.txt');



pos0 = [40.0678531*glv.deg; 116.247326*glv.deg; 47.92];     %��ʼλ��   
vn0 = [0; 0; 0];                                            %��ʼ�ٶȣ���ʱΪ��ֹ̬
att0 = [0.0040; -0.0048; -2.9976];                          %����UKF�˲���׼���׶ε���ֵ̬
cnb0 = a2cnb(att0);
qnb0 = m2qnb(cnb0);                                         %��ֵ̬�仯Ϊ��ʼ��Ԫ��
qnb = qnb0; vn = vn0; pos = pos0;
n = 1;                                                      %n������
Ts = n*Ts;                                                  %��������
dpos = [1./glv.Re; 1./glv.Re; 0.3];                         %λ�ù۲�����
dv = [.5; .5; .1];                                          %�ٶȹ۲�����
[eb, web, db, wdb] = drift_bias([1;1;1;1]);                 %����������ƫ�����Ư��

%%%%%%%%%%%%%�˲����������á���ѡȡ����Ϊ�˲�����֪ʶ��ʵ�ʹ�������%%%%%%%%%%%%%%%
Qt = diag([web; wdb; repmat(0.01, 3, 1); 0.1*web; 0.1*wdb])^2;
Rk = diag([dv; dpos])^2;
Pk = diag([[.5;0.5;2.]*glv.deg; [2.; 2.; 1.]; ...
    [5./glv.Re; 5./glv.Re; 2.]; [2.; 2.; 2.]*glv.dph; ...
    [3.; 3.; 3.]*glv.mg])^2;
Xk = [[0.05;0.05;0.5]*glv.deg; [.3; .3; .1];...
    [.5/glv.Re; .5/glv.Re; .1]; [.7; .7; .7]*glv.dph;...
    [1.; 1.; 1.]*glv.mg];
Hk = [zeros(6, 3), eye(6), zeros(6,6)];
len = 262250;                                               %����ѭ����������
lenp = 1320;
Xkk = zeros(lenp, length(Xk));                              %�˲�״̬�洢
Zk_kalman = zeros(lenp, 6);                                 %�۲����洢�����ٶȺ�λ�����
navigation = zeros(len, 9);                                 %���������洢
kv = 2035;                                                  %��˼����е�IMU���ݶ���
kk = 1;
kerr = 1;



for k = 454321:n:716560                               %10.5�����г����룬20����

   k1 = k+n-1;                                             %�г�ʱ��
   wm = imu0(k:k1,6:8)';
   vm = imu0(k:k1,3:5)';
   [qnb, vn, pos] = sinss(qnb, vn, pos, wm, vm, Ts);
   cnb = q2cnb(qnb);
   [x, y, z] = DCMtoEuler(cnb);
   att = [x; y; z];
   navigation(kk,:) = [att'./glv.deg, vn', pos(1:2,1)'/glv.deg, pos(3,1)];
   wm = sum(wm,2);                                          %����õ���������Ӧ
   vm = sum(vm,2);                                          %getfsҪ����������ڵ�����
   Ft = getfs(qnb, vn, pos, vm./Ts, 15);                    %����ϵͳϵ������
   [Fikk_1, Qk] = kfdis(Ft, Qt, Ts, 4);                     %״̬һ��ת�ƾ���
   for j = k:1:k1
       if imu0(j, 2)==veli(kv, 2)                           %��GPS��������в���        
           vi = veli(kv, 3:5)';
           pi = [posi(kv, 3:4)'*glv.deg; posi(kv, 5)];
           Zk = [vn-vi; pos-pi];                            %����۲���
           Zk_kalman(kerr, :) = Zk';
           [Xk, Pk, Kk] = kalman(Fikk_1, Qk, Xk, Pk, Hk, Rk, Zk);
           Xkk(kerr, :) = Xk';
           vn = vn-Xk(4:6, 1);                              %�γɱջ�����
           pos = pos-Xk(7:9, 1);                            %�ջ�����
           pos(3) = posi(kv, 5);                            %���Ƹ߶�ͨ����ɢ
           xk(9)=0;                                         %���Ƹ߶�ͨ����ɢ
           kv = kv+1;
           kerr = kerr+1;
       else
           [Xk, Pk] = kalman(Fikk_1, Qk, Xk, Pk);           %�޹۲���ʱֻ����Ԥ��
       end
   end
   timedisp(kk, Ts, 100);                                   %��ʾ����
   kk = kk+1;
end

save navigation_kalman navigation
save Xkk_kalman Xkk
save Zk_kalman Zk_kalman

 figure,             %������ϵ����켣
 plot(navigation(1:kk-1,8), navigation(1:kk-1,7), '-'), grid on
xlabel('\it\lambda\rm / \circ'); ylabel('\itLatitude\rm / \circ'); 

figure,             %����������������
plot(Xkk(1:kerr-1,1)/glv.deg, '-b'), grid on
xlabel('\itt/s'); ylabel('\it\theta\rm / \circ'); grid on

figure,             %����������������
plot(Xkk(1:kerr-1,2)/glv.deg, '-b'), grid on
xlabel('\itt/s'); ylabel('\it\gamma\rm / \circ');

figure,             %��������λ�������
plot(Xkk(1:kerr-1,3)/glv.deg, '-b'), grid on
xlabel('\itt/s'); ylabel('\it\psi\rm / \circ');

figure,             %�����������ٶ����
subplot(2,1,1), plot(Zk_kalman(1:kerr-1,1)), grid on
xlabel('\itt/s'); ylabel('\it\delta Ve\rm /m/s'); 
subplot(2,1,2), plot(Zk_kalman(1:kerr-1,2)), grid on
xlabel('\itt/s'); ylabel('\it\delta Vn\rm /m/s'); 

figure,             %����γ�ȡ��������
subplot(2,1,1), plot(Zk_kalman(1:kerr-1,4)*glv.Re), grid on
xlabel('\itt/s'); ylabel('\it\delta Latitude\rm /m'); 
subplot(2,1,2), plot(Zk_kalman(1:kerr-1,5)*glv.Re), grid on
xlabel('\itt/s'); ylabel('\it\delta lambda\rm /m'); 

toc

