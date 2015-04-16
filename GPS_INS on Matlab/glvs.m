% function glvs
%全局变量
    global glv
    glv.Re = 6378137;               %地球半径(GPS-84)
    glv.f = 1/298.257;              %地球扁率
    glv.Rp = (1-glv.f)*glv.Re;      %地球椭圆度等其它几何参数
    glv.e = sqrt(2*glv.f-glv.f^2); glv.e2 = glv.e^2;%偏心率
    glv.ep = sqrt(glv.Re^2+glv.Rp^2)/glv.Rp;	glv.ep2 = glv.ep^2;
    glv.wie = 7.2921151467e-5;      %地球自转角速率
    glv.meru = glv.wie/1000;        %毫地转率
    glv.g0 = 9.7803267714;          %重力加速度
    glv.mg = 1.0e-3*glv.g0;         %毫重力加速度
    glv.ug = 1.0e-6*glv.g0;         %微重力加速度
    glv.ugpg2 = glv.ug/glv.g0^2;    %加表二次项系数
    glv.ppm = 1.0e-6;               %百万分之一
    glv.deg = pi/180;               %角度
    glv.min = glv.deg/60;           %角分
    glv.sec = glv.min/60;           %角秒
    glv.hur = 3600;                 %小时
    glv.dps = pi/180/1;             %度/秒
    glv.dph = glv.deg/glv.hur;      %度/小时
    glv.dpss = pi/180/sqrt(1);          %度/sqrt(秒)
    glv.dpsh = glv.deg/sqrt(glv.hur);   %度/sqrt(小时), 角度随机游走系数
    glv.Hz = 1/1;                   %赫兹(1/s)
    glv.ugpsHz = glv.ug/sqrt(glv.Hz);    %ug/sqrt(Hz), 速度随机游走系数
    glv.mil = 2*pi/6000;            %密位
    glv.nm = 1853;                  %海里
    glv.kn = glv.nm/glv.hur;        %节
    glv.cs = [                      %圆锥和划船误差补偿系数
        2./3,       0,          0,          0
        9./20,      27./20,     0,          0
        54./105,    92./105,    214./105,    0
        250./504,   525./504,   650./504,   1375./504 ];

