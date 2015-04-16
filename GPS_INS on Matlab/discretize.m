function g=discretize(f,t)
%g=discretize(f,t)
% f:连续函数  t:采样间隔
k1=@(x)f(x);
k2=@(x)f(x+k1(x)*t/2);
k3=@(x)f(x+k2(x)*t/2);
k4=@(x)f(x+k3(x)*t);
g=@(x)(k1(x)+2*k2(x)+2*k3(x)+k4(x))*t/6+x;%四阶龙格库塔法