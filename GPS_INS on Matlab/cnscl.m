function [phim, dvbm] = cnscl(wm, vm)
global glv
    n = size(wm,2);                                   %n为列数
    cm = [0;0;0]; sm = [0;0;0];
    for k=1:n-1
        cm = cm + glv.cs(n-1,k)*wm(:,k);              %准备圆锥误差补偿
    end
    wmm = sum(wm,2);                                  %行和的列向量
    phim = wmm + cross(cm,wm(:,n));
    if nargin==2
        for k=1:n-1
            sm = sm + glv.cs(n-1,k)*vm(:,k);          %准备划船误差补偿
        end
        vmm = sum(vm,2);
	    dvbm = vmm + 1.0/2*cross(wmm,vmm) + (cross(cm,vm(:,n))+cross(sm,wm(:,n)));
    end
