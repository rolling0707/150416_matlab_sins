function att=m2att(Cnb)                    %×ËÌ¬¾ØÕóÇóĞı×ª½Ç
    att = [ asin(Cnb(3,2));                %¸©Ñö½Ç
            atan2(-Cnb(3,1),Cnb(3,3));     %ºá¹ö½Ç
            atan2(-Cnb(1,2),Cnb(2,2)) ];   %º½Ïò½Ç
