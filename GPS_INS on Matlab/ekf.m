function [Xk,Pk,Kk] = ekf(Phikk_1,Xk_1,Qk,fXk_1,Pk_1,Hk,Rk,Zk, Ts)
%%%参考《卡尔曼滤波与组合导航原理》p186,6.2.88-6.2.93    
    Xkk_1=Xk_1+fXk_1*Ts;
    Pkk_1=Phikk_1*Pk_1*Phikk_1'+Qk;
    Kk=Pkk_1*Hk'*(Hk*Pkk_1*Hk'+Rk)^-1;
%     Pk=Pkk_1-Kk*(Hk*Pkk_1*Hk'+Rk)*Kk';
    Pk=(eye(15)-Kk*Hk)*Pkk_1*(eye(15)-Kk*Hk)'+Kk*Rk*Kk';
    detalXk=Kk*(Zk-Hk*Xkk_1);
    Xk=Xkk_1+detalXk;
%     Pkk_1=Phikk_1*Pk_1*Phikk_1'+Qk;
%     Kk=Pkk_1*Hk'*(Hk*Pkk_1*Hk'+Rk)^-1;
%     Xk=fXk_1+Kk*Zk_hfX;
%     Pk=Pkk_1-Kk*(Hk*Pkk_1*Hk'+Rk)*Kk';