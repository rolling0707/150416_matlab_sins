function [Xk, Pk, Kk] = kalman(Phikk_1, Qk, Xk_1, Pk_1, Hk, Rk, Zk)
    

if nargin<7     % 无量测时仅进行状态递推(时间更新)
        
	Xk = Phikk_1*Xk_1;
        
	Pk = Phikk_1*Pk_1*Phikk_1'+Qk;
    
else            % 有测量时滤波
        
	Xkk_1 = Phikk_1*Xk_1;    
        
	Pkk_1 = Phikk_1*Pk_1*Phikk_1' + Qk; 
        
	Pxz = Pkk_1*Hk';    
	Pzz = Hk*Pxz + Rk;     
	Kk = Pxz*Pzz^-1;
        
	Xk = Xkk_1 + Kk*(Zk-Hk*Xkk_1);

	%Pk = Pkk_1 - (Kk*Pxz')'; 
	%Pk = (Pk+Pk')/2; 

	%Pk = Pkk_1 - Pxz*Kk'; 
	%Pk = (Pk+Pk')/2; 
        
	Pk = Pkk_1 - Kk*Pzz*Kk';
    
end
