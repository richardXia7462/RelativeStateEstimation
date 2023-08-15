function DrawEstimation3Sigma(t, error, cov, bPVQ, bBiasF, bBiasL)

    if bPVQ
        errAtt = error{1}; errTra = error{2}; errVel = error{3};
        covAtt = cov{1}; covTra = cov{2}; covVel = cov{3};
        tit='attitude errors'; ylab={'x (rad)';'y (rad)';'z (rad)'};
        PlotEstSig(t,errAtt,covAtt,ylab,tit,'r','b:');
        tit='position errors'; ylab={'x (m)';'y (m)';'z (m)'};
        PlotEstSig(t,errTra,covTra,ylab,tit,'r','b:')
        tit='velocity errors'; ylab={'x (m/s)';'y (m/s)';'z (m/s)'};
        PlotEstSig(t,errVel,covVel,ylab,tit,'r','b:')
    end
    
    if bBiasF
        errBiasGF = error{4}; errBiasAF = error{5}; 
        covBiasGF = cov{4}; covBiasAF = cov{5}; 
        tit='gyro bias errors'; ylab={'x (rad/s)';'y (rad/s)';'z (rad/s)'};
        PlotEstSig(t,errBiasGF,covBiasGF,ylab,tit,'r','b:')
        tit='acc bias errors'; ylab={'x (m/s^2)';'y (m/s^2)';'z (m/s^2)'};
        PlotEstSig(t,errBiasAF,covBiasAF,ylab,tit,'r','b:')
    end
    
    if bBiasL
        errBiasGL = error{6}; errBiasAL = error{7}; 
        covBiasGL = cov{6}; covBiasAL = cov{7}; 
        tit='leader gyro bias errors'; ylab={'x (rad/s)';'y (rad/s)';'z (rad/s)'};
        PlotEstSig(t,errBiasGL,covBiasGL,ylab,tit,'r','b:')
        tit='leader acc bias errors'; ylab={'x (m/s^2)';'y (m/s^2)';'z (m/s^2)'};
        PlotEstSig(t,errBiasAL,covBiasAL,ylab,tit,'r','b:')
    end






        

