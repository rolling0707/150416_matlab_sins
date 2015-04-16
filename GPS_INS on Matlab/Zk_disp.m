clear
glvs
load('Zk_ekf.mat');
load('Zk_ukf.mat');
figure,
subplot(2,1,1), plot(Zk_ekf(1:1262, 1), '-b'); grid on
hold on
plot(Zk_ukf(1:1262, 1), '-r');

subplot(2,1,2), plot(Zk_ekf(1:1262, 2), '-b'); grid on
hold on
plot(Zk_ukf(1:1262, 2), '-r');

figure,
subplot(2,1,1), plot(Zk_ekf(1:1262, 4)*glv.Re, '-b'); grid on
hold on
plot(Zk_ukf(1:1262, 4)*glv.Re, '-r');

subplot(2,1,2), plot(Zk_ekf(1:1262, 5)*glv.Re, '-b'); grid on
hold on
plot(Zk_ukf(1:1262, 5)*glv.Re, '-r');