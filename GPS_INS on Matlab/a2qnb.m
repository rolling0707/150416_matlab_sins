function qnb = a2qnb(att) %由旋转角求四元数
    qnb = m2qnb(a2cnb(att));