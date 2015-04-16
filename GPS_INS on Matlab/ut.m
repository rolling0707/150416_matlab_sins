function [y,Y,P,Y1] = ut(f,X,Wm,Wc,n,R)
    L=size(X,2);%sigma点的个数--》2n+1
    y=zeros(n,1);%系统变量或者观测变量的初始化
    Y=zeros(n,L);%sigma点一步预测的初始化
    for k=1:L
%         Y(:,k)=f(X(:,k),kx,ky,g,Ts);%sigma点或者观测量的一步预测
        Y(:,k) = f(X(:,k));
        y=y+Wm(k)*Y(:,k);%sigma点或者观测量的求平均
    end
    Y1=Y-y(:,ones(1,L));%一步预测和真实值的误差向量
    P=Y1*diag(Wc)*Y1'+R;%系统变量协方差的一步预测


    