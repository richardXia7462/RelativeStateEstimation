function dst=GetCrossMat(src)
    dst=[0 -src(3) src(2);
        src(3) 0 -src(1);
        -src(2) src(1) 0];
end