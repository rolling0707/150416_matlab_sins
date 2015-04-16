function [eb, web, db, wdb] = drift_bias(val)
% 设置陀螺和加计的误差：常值漂移与随机游走系数
% 输入 --- 一般为4X3矩阵，第一行为陀螺常值漂移，乘单位0.01deg/hur
%                         第二行为陀螺角度游走系数，乘单位0.001deg/sqrt(hur)
%                         第三行为加计常值漂移，乘单位100ug
%                         第四行为加计速度游走系数，乘单位10ug/sqrt(Hz)
%           如果为4X1向量，则每种类型误差的各轴数值相同，再乘以上述单位
%           如果为标量，则所有误差数值相同，再乘以上述单位
% 输出 ---
%     eb - 陀螺常值漂移    web - 陀螺随机游走系数
%     db - 加计常值漂移    wdb - 加计随机游走系数
    global glv
    if exist('val')==0
        val = 1;
    end
    [m, n] = size(val);
    if m == 1
        val = val*ones(4,1);
    end
    if n == 1
        val = repmat(val,1,3);
    end
    eb = val(1,:)'*0.5*glv.dph;
    web = val(2,:)'*0.08*glv.dpsh;
    db = val(3,:)'*800*glv.ug;
    wdb = val(4,:)'*35*glv.ugpsHz;