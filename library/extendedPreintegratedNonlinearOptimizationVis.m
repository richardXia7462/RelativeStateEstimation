function [x_est, p_cov] = extendedPreintegratedNonlinearOptimizationVis(IMUF, IMUL, imuParameter, imuParameterL, F, Tcb, camParameters, bPixel, obsPixelPH, rPixel, x0, p0)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

nS = size(obsPixelPH,1);
CovGyro = imuParameterL(1)^2*eye(6);
x_est = zeros(nS,16);
p_cov = zeros(nS,15);

tPre = obsPixelPH(1,1);
[imuF, iteIMUF, imuL, iteIMUL] = GetIMUAtFirstObservation(tPre, IMUF, IMUL);

Rlf0 = quat2rotm(x0(1:4)); 
tlf0 = x0(5:7)';
vlf0 = x0(8:10)'; 
bg0 = x0(11:13)'; 
ba0 = x0(14:16)';
infoPrior0 = inv(p0(1:15,1:15));
x_est(1,:) = x0(1:16);
p_cov(1,:) = diag(p0(1:15,1:15));
nFace = size(bPixel,2);
nTagPerFace = 4;
for i = 2:nS
    % i = i+1
    tCur = obsPixelPH(i,1);
    
    j = iteIMUF;
    while j<size(IMUF,1) && IMUF(j,1)<tCur
        j = j+1;
    end
    imuF = [imuF; IMUF(iteIMUF:j-1,:)];
    tab = IMUF(j,1) - IMUF(j-1,1);
    tini = tCur - IMUF(j-1,1);
    angVel = IMUF(j-1,2:4) + (IMUF(j,2:4) - IMUF(j-1,2:4)) * (tini / tab);
    acc = IMUF(j-1,5:7) + (IMUF(j,5:7) - IMUF(j-1,5:7)) * (tini / tab);
    imuF = [imuF; tCur, angVel, acc];
    iteIMUF = j;

    j = iteIMUL;
    while j<size(IMUL,1) && IMUL(j,1)<tCur
        j = j+1;
    end
    imuL = [imuL; IMUL(iteIMUL:j-1,:)];

    tab = IMUL(j,1) - IMUL(j-1,1);
    tini = tCur - IMUL(j-1,1);
    angVel = IMUL(j-1,2:4) + (IMUL(j,2:4) - IMUL(j-1,2:4)) * (tini / tab);
    acc = IMUL(j-1,5:7) + (IMUL(j,5:7) - IMUL(j-1,5:7)) * (tini / tab);
    imuL = [imuL; tCur, angVel, acc];
    iteIMUL = j;

    [del_a_f,del_v_f,del_p_f,CovF,JacoInt,CovG,CovA] = preIntegrateMeasurement2(imuF(:,5:7),...
        imuF(:,2:4),imuF(:,1),imuParameter,bg0,ba0);
    [del_a_l,del_v_l,del_p_l,CovL,~,~,~,correl] = preIntegrateMeasurement2(imuL(:,5:7),...
        imuL(:,2:4),imuL(:,1),imuParameter,zeros(3,1),zeros(3,1)); 
    gyrol0 = imuL(1,2:4)'; gyrol1 = imuL(end,2:4)';       
    imuF = imuF(end,:);
    imuL = imuL(end,:);     
    T = tCur-tPre;
    tPre = tCur;

    tlf = del_a_l'*(Rlf0*del_p_f-del_p_l+tlf0+(vlf0+x_crossMat(gyrol0)*tlf0)*T);
    vlf = del_a_l'*(Rlf0*del_v_f-del_v_l+vlf0+x_crossMat(gyrol0)*tlf0)-x_crossMat(gyrol1)*tlf;
    Rlf = del_a_l'*Rlf0*del_a_f;  
    bg = bg0; ba = ba0;

    Mr = del_a_l'*Rlf0*del_a_f;
    AF = [eye(3), zeros(3), zeros(3);
        zeros(3), Rlf0, zeros(3);
        zeros(3), zeros(3), Rlf0];
    AL = [-Mr', zeros(3), zeros(3);
        del_a_l*x_crossMat(vlf+x_crossMat(gyrol1)*tlf), -eye(3), zeros(3);
        del_a_l*x_crossMat(tlf), zeros(3), -eye(3)];
    BL = [zeros(3), zeros(3);
        -x_crossMat(tlf0), del_a_l*x_crossMat(tlf);
        -x_crossMat(tlf0)*T, zeros(3)];
    S = AL*correl*BL(:,1:3)'; S = S+S';
    infoPreIntegration = inv(AF*CovF*AF'+AL*CovL*AL'+BL*CovGyro*BL'+S);
    infoBiasGyro = inv(CovG); infoBiasAcc = inv(CovA);

    r_est = []; h = []; obs=[];
    for j = 1:nFace
        if bPixel(i, j)==1
            obs = [obs; obsPixelPH(i,(j-1)*nTagPerFace*2+2:j*nTagPerFace*2+1)'];
        end
    end
    infoVis = 1/rPixel*eye(length(obs));

    a_f2l0_mea = Rlf0; p_f2l0_mea = tlf0; v_f2l0_mea = vlf0; bg0_mea = bg0; ba0_mea = ba0;    

    infoAll = zeros(30+length(obs),30+length(obs));
    infoAll(1:9,1:9) = infoPreIntegration;
    infoAll(10:24,10:24) = infoPrior0;            
    infoAll(25:27,25:27) = infoBiasGyro;
    infoAll(28:30,28:30) = infoBiasAcc;
    infoAll(31:end, 31:end) = infoVis;
    
    Norm = 0;
    for ite=1:1  
        Mr = del_a_l'*Rlf0*del_a_f;
        dMr = x_LieLog(Rlf'*Mr)';
        dMv = Rlf0*del_v_f-del_v_l+vlf0+x_crossMat(gyrol0)*tlf0-del_a_l*(vlf+x_crossMat(gyrol1)*tlf);
        dMp = Rlf0*del_p_f-del_p_l-del_a_l*tlf+tlf0+(vlf0+x_crossMat(gyrol0)*tlf0)*T;
        errPreIntegration = [dMr; dMv; dMp];
        jPreInt = preIntegrateMeasurementErrJacobian(del_a_f,del_v_f,del_p_f,JacoInt,Rlf0,Rlf,del_a_l,gyrol0,gyrol1,T);      

        r_est = []; h = []; 
        for j = 1:nFace
            if bPixel(i, j)==1
                rEstTemp = zeros(2*nTagPerFace,1);
                hTemp = zeros(2*nTagPerFace,30);
                for k = 1:nTagPerFace
                    v3D = Tcb(1:3,:)*[Rlf*F(:,(j-1)*nTagPerFace+k)+tlf; 1];
                    rEstTemp((k-1)*2+1:k*2) = Project_PinHole(v3D, camParameters);
                    Jac = projectJac_PinHole(v3D, camParameters);             
                    hTemp((k-1)*2+1:k*2,1:3) = -Jac*Tcb(1:3,1:3)*Rlf*x_crossMat(F(:,(j-1)*nTagPerFace+k));
                    hTemp((k-1)*2+1:k*2,4:6) = Jac*Tcb(1:3,1:3);
                end
                r_est = [r_est; rEstTemp];
                h = [h; hTemp];
            end
        end
        errFeature = obs-r_est;
        jacFeature = -h;

        priorState = [rotm2quat(a_f2l0_mea)';p_f2l0_mea;v_f2l0_mea;bg0_mea;ba0_mea];
        estState = [rotm2quat(Rlf0)';tlf0;vlf0;bg0;ba0];
        [errPrior,jPrior] = priorErrJac(priorState,estState);

        [errBiasGyro,jGyro,errBiasAcc,jAcc] = biasErrJac(bg,bg0,ba,ba0);

        Jacobian = [jPreInt; jPrior; jGyro; jAcc; jacFeature];
        errAll = [errPreIntegration; errPrior; errBiasGyro; errBiasAcc; errFeature];

        
        if Norm~=0 && abs(norm(errAll)-Norm)/Norm<2*1e-3
            break
        end
%         Norm = norm(errAll)
        deltaX = -(Jacobian'*infoAll*Jacobian)\Jacobian'*infoAll*errAll; %norm(deltaX) 

        dqj = deltaX(1:3); dtj = deltaX(4:6); dvj = deltaX(7:9);        
        dbgj = deltaX(10:12); dbaj = deltaX(13:15); 
        dqi = deltaX(16:18); dti = deltaX(19:21); dvi = deltaX(22:24);       
        dbgi = deltaX(25:27); dbai = deltaX(28:30); 

        Rlf = Rlf*x_LieExp(dqj); tlf = tlf+dtj; vlf = vlf+dvj; 
        Rlf0 = Rlf0*x_LieExp(dqi); tlf0 = tlf0+dti; vlf0 = vlf0+dvi;
        bg = bg+dbgj; ba = ba+dbaj;
        bg0 = bg0+dbgi; ba0 = ba0+dbai;
        [del_a_f,del_v_f,del_p_f] = getDeltaIntgrateMeasurement(del_a_f,del_v_f,del_p_f,JacoInt,dbgi,dbai);
    end

    Rlf0 = Rlf; tlf0 = tlf; vlf0 = vlf; bg0 = bg; ba0 = ba;       
    Info = Jacobian'*infoAll*Jacobian;
    H11 = Info(1:15,1:15); H10 = Info(1:15,16:30); H01 = Info(16:30,1:15); H00 = Info(16:30,16:30);
    infoPrior0 = H11-H10*inv(H00)*H01;
    x_est(i,:) = [rotm2quat(Rlf), tlf', vlf', bg', ba'];
    p_cov(i,:) = diag(inv(infoPrior0));
    
end
end

