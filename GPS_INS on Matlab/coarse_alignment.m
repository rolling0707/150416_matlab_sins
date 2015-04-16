function [pitch, roll, yaw] = coarse_alignment(data, pos)
%%%算法参考《捷联惯导系统粗对准方法比较》-魏春玲，张洪越-2000
    glvs
    Ts = 0.005;
%     pos = [40.067987*glv.deg; 116.247871*glv.deg; 43.53];
    [wnie, wnen, gn, retp] = earth(pos);
    %Lat = 40.0678524*pi/180;
    len = length(data);
    %wnie = glv.wie*[0; cos(Lat); sin(Lat)];
    gn = [0; 0; -retp.g];
    rn = askew(gn)*wnie;
    rrn = askew(rn)*gn;
    
    wbie_measum = [0; 0; 0];
    gb_measum = [0; 0; 0];
    for k=1:1:len
        [phim, dvbm] = cnscl(data(k, 6:8)', data(k, 3:5)');
        wbie_measum = wbie_measum+phim;
        gb_measum = gb_measum+dvbm;
    end
    wbie_mean = wbie_measum/len/Ts;
    gb_mean = -gb_measum/len/Ts;
    rb = askew(gb_mean)*wbie_mean;%-gb_mean为gn在机体系下的分量
    rrb = askew(rb)*gb_mean;
    
    cnb_align = ([gn'; rn'; rrn'])\([gb_mean'; rb'; rrb']);%（公式8）
    Cnb = cnb_align/(sqrtm(cnb_align'*cnb_align));%（公式7）正交化处理
    [pitch, roll, yaw] = DCMtoEuler(Cnb);
  