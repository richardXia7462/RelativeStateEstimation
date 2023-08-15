function [mpImuFromLastFrame] = ImuSynchronization(imuL,imuF,tCur)
% imuL = vImuMeasL; imuF = vImuMeas;
%UNTITLED2 此处显示有关此函数的摘要
%   mpImuFromLastFrame = [t gyroF accF gyroL accL]

    nF = size(imuF,1);
    nL = size(imuL,1);
    if nL<2 || nF<2
        error('imu not enough')
    end
%     if imuF(1,1)>tPre || imuL(1,1)>tPre
%         error('first imu wrong')
%     end
    if imuF(nF,1)<tCur || imuL(nL,1)<tCur
        warning('last imu wrong')
    end
    if imuF(nF-1,1)>tCur || imuL(nL-1,1)>tCur
        error('last but one imu wrong')
    end
    
    mpImuFromLastFrame = [];    
    ite = 1;
    for i=1:nL-1
        % i = i+1;
        if imuL(i,1)>=imuF(ite,1)
            while ite+1<=nF && imuL(i,1)>=imuF(ite+1,1)
                ite = ite+1;
            end
            if ite>=nF
                break
            end
            t = imuF(ite,1); t_ = imuF(ite+1,1);
            tab = t_ - t;
            tini = imuL(i,1) - t;
            angVel = imuF(ite,2:4) + (imuF(ite+1,2:4) - imuF(ite,2:4)) * (tini / tab);
            acc = imuF(ite,5:7) + (imuF(ite+1,5:7) - imuF(ite,5:7)) * (tini / tab);
            mpImuFromLastFrame = [mpImuFromLastFrame; imuL(i,:), angVel, acc];
        end
    end
    
	t = imuL(nL-1,1); t_ = imuL(nL,1);   
    tab = t_ - t;
    tini = tCur - t;
    angVel = imuL(nL-1,2:4) + (imuL(nL,2:4) - imuL(nL-1,2:4)) * (tini / tab);
    acc = imuL(nL-1,5:7) + (imuL(nL,5:7) - imuL(nL-1,5:7)) * (tini / tab);
    imu = [tCur, angVel, acc];

	t = imuF(nF-1,1); t_ = imuF(nF,1);
    tab = t_ - t;
    tini = tCur - t;
    angVel = imuF(nF-1,2:4) + (imuF(nF,2:4) - imuF(nF-1,2:4)) * (tini / tab);
    acc = imuF(nF-1,5:7) + (imuF(nF,5:7) - imuF(nF-1,5:7)) * (tini / tab);
    mpImuFromLastFrame = [mpImuFromLastFrame; imu, angVel, acc];
end

