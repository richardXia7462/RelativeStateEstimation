function [x_est, p_cov] = NLOpt_GN_PH_ExtPreint_BothBias(IMUF, IMUL, x0, p0, obsPixelPH, rPixel, bPixel, markF, Tcb, camParameters, imuParameter, imuParameterL)
 
nFace = size(bPixel,2);
nTagPerFace = (size(obsPixelPH,2)-1)/nFace/2;
num = size(obsPixelPH,1);
x_est = zeros(num,22);
p_cov = zeros(num,21);
x_est(1,:) = x0;
p_cov(1,:) = diag(p0);
infoPrior0 = inv(p0);
CovGyro = imuParameterL(1)^2*eye(6);

tPre = obsPixelPH(1,1);
[imuF, iteIMUF, imuL, iteIMUL] = GetIMUAtFirstObservation(tPre, IMUF, IMUL);

Rlf0 = quat2rotm(x0(1:4)); 
tlf0 = x0(5:7)';
vlf0 = x0(8:10)'; 
bgf0 = x0(11:13)'; 
baf0 = x0(14:16)';
bgl0 = x0(17:19)'; 
bal0 = x0(20:22)';
for i = 2:num

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
        
    [del_a_f,del_v_f,del_p_f,CovF,JacoInt,CovG,CovA] = PreIntegrateMeasurementBiases(imuF(:,5:7),...
        imuF(:,2:4),imuF(:,1),imuParameter,bgf0,baf0);
    [del_a_l,del_v_l,del_p_l,CovL,JacoIntL,CovGL,CovAL,correl] = PreIntegrateMeasurementBiases(imuL(:,5:7),...
        imuL(:,2:4),imuL(:,1),imuParameterL,bgl0,bal0); 
    gyrol0 = imuL(1,2:4)'-bgl0; gyrol1 = imuL(end,2:4)'-bgl0;       
    imuF = imuF(end,:);
    imuL = imuL(end,:);     
    tij = tCur-tPre;
    tPre = tCur;
    
    % propogate state
    tlf = del_a_l'*(Rlf0*del_p_f-del_p_l+tlf0+(vlf0+GetCrossMat(gyrol0)*tlf0)*tij);
    vlf = del_a_l'*(Rlf0*del_v_f-del_v_l+vlf0+GetCrossMat(gyrol0)*tlf0)-GetCrossMat(gyrol1)*tlf;
    Rlf = del_a_l'*Rlf0*del_a_f;  
    bgf = bgf0; baf = baf0;
    bgl = bgl0; bal = bal0;
    
    AF = [eye(3), zeros(3), zeros(3);
        zeros(3), Rlf0, zeros(3);
        zeros(3), zeros(3), Rlf0];
    Mr = del_a_l'*Rlf0*del_a_f;
    AL = [-Mr', zeros(3), zeros(3);
        del_a_l*GetCrossMat(vlf+GetCrossMat(gyrol1)*tlf), -eye(3), zeros(3);
        del_a_l*GetCrossMat(tlf), zeros(3), -eye(3)];
    BL = [zeros(3), zeros(3);
        -GetCrossMat(tlf0), del_a_l*GetCrossMat(tlf);
        -GetCrossMat(tlf0)*tij, zeros(3)];
%     infoPreIntegration = inv(AF*CovF*AF'+AL*CovL*AL'+BL*CovGyro*BL');
    S = AL*correl*BL(:,1:3)'; S = S+S';
    infoPreIntegration = inv(AF*CovF*AF'+AL*CovL*AL'+BL*CovGyro*BL'+S);
    infoBiasGyro = inv(CovG); infoBiasAcc = inv(CovA);
    infoBiasGyroL = inv(CovGL); infoBiasAccL = inv(CovAL);
    
    Rlf0_mea = Rlf0; tlf0_mea = tlf0; 
    vlf0_mea = vlf0; 
    bgf0_mea = bgf0; baf0_mea = baf0;    
    bgl0_mea = bgl0; bal0_mea = bal0;
    
    nFeaDim = sum(bPixel(i,:))*4*2;
    infoAll = eye(42+nFeaDim,42+nFeaDim);
    index = 1;     
    infoAll(index:index+8,index:index+8) = infoPreIntegration; 
    index = index+9;
    infoAll(index:index+20,index:index+20) = infoPrior0;  
    index = index+21;
    infoAll(index:index+2,index:index+2) = infoBiasGyro;  
    index = index+3;
    infoAll(index:index+2,index:index+2) = infoBiasAcc;  
    index = index+3;
    infoAll(index:index+2,index:index+2) = infoBiasGyroL;
    index = index+3;
    infoAll(index:index+2,index:index+2) = infoBiasAccL;
    index = index+3;
    infoAll(index:index+nFeaDim-1,index:index+nFeaDim-1) = 1/rPixel*eye(nFeaDim);
    
    for ite = 1:1
        Mr = del_a_l'*Rlf0*del_a_f;
        dMr = LieLog(Rlf'*Mr)';
        dMv = Rlf0*del_v_f-del_v_l+vlf0+GetCrossMat(gyrol0)*tlf0-del_a_l*(vlf+GetCrossMat(gyrol1)*tlf);
        dMp = Rlf0*del_p_f-del_p_l-del_a_l*tlf+tlf0+(vlf0+GetCrossMat(gyrol0)*tlf0)*tij;
        errPreIntegration = [dMr; dMv; dMp];    
        jPreInt = PreIntegrateMeasurementErrJacobianBias(del_a_f,del_v_f,del_p_f,JacoInt,JacoIntL,Rlf0,Rlf,vlf,tlf0,tlf,del_a_l,gyrol0,gyrol1,tij);
        
        priorState = [rotm2quat(Rlf0_mea)'; tlf0_mea; vlf0_mea; bgf0_mea; baf0_mea; bgl0_mea; bal0_mea];
        estState = [rotm2quat(Rlf0)'; tlf0; vlf0; bgf0; baf0; bgl0; bal0];
        [errPrior,jPrior] = PriorErrJacBias(priorState,estState);
        
        [errBiasGyro,jGyro_,errBiasAcc,jAcc_] = BiasErrJacBias(bgf,bgf0,baf,baf0, tij);
        jacGyro = [zeros(3,9) jGyro_{1} zeros(3,18) jGyro_{2} zeros(3,9)];
        jacAcc = [zeros(3,12) jAcc_{1} zeros(3,18) jAcc_{2} zeros(3,6)];
        [errBiasGyroL,jGyroL,errBiasAccL,jAccL] = BiasErrJacBias(bgl, bgl0, bal, bal0, tij);
        jacGyroL = [zeros(3,15) jGyroL{1} zeros(3,3) zeros(3,15) jGyroL{2} zeros(3,3)];
        jacAccL = [zeros(3,18) jAccL{1} zeros(3,18) jAccL{2}];
        
        % feature observation
        obs=[];
        for j = 1:nFace
            if bPixel(i, j)==1
                obs = [obs; obsPixelPH(i,(j-1)*nTagPerFace*2+2:j*nTagPerFace*2+1)'];
            end
        end
        rcov = rPixel*eye(length(obs));
        
        r_est = []; h = []; 
        for j = 1:nFace
            if bPixel(i, j)==1
                rEstTemp = zeros(2*nTagPerFace,1);
                hTemp = zeros(2*nTagPerFace,42);
                for k = 1:nTagPerFace
                    v3D = Tcb(1:3,:)*[Rlf*markF(:,(j-1)*nTagPerFace+k)+tlf; 1];
                    rEstTemp((k-1)*2+1:k*2) = ProjectPinHole(v3D, camParameters);
                    Jac = ProjectJacPinHole(v3D, camParameters);             
                    hTemp((k-1)*2+1:k*2,1:3) = -Jac*Tcb(1:3,1:3)*Rlf*GetCrossMat(markF(:,(j-1)*nTagPerFace+k));
                    hTemp((k-1)*2+1:k*2,4:6) = Jac*Tcb(1:3,1:3);
                end
                r_est = [r_est; rEstTemp];
                h = [h; hTemp];
            end
        end
        errFeature = obs-r_est;
        jacFeature = -h;
        
        Jacobian = [jPreInt; jPrior; jacGyro; jacAcc; jacGyroL; jacAccL; jacFeature];
        errAll = [errPreIntegration; errPrior; errBiasGyro; errBiasAcc; errBiasGyroL; errBiasAccL; errFeature];
        deltaX = -(Jacobian'*infoAll*Jacobian)\Jacobian'*infoAll*errAll;

        dqj = deltaX(1:3); dtj = deltaX(4:6); dvj = deltaX(7:9); 
        dgfj = deltaX(10:12); dafj = deltaX(13:15);
        dbglj = deltaX(16:18); dbalj = deltaX(19:21); 
        indState = 21;

        dqi = deltaX((1:3)+indState); dti = deltaX((4:6)+indState); dvi = deltaX((7:9)+indState);    
        dbgfi = deltaX((10:12)+indState); dbafi = deltaX((13:15)+indState);
        dbgli = deltaX((16:18)+indState); dbali = deltaX((19:21)+indState); 
        Rlf = Rlf*LieExp(dqj); tlf = tlf+dtj; vlf = vlf+dvj; 
        bgf = bgf+dgfj; baf = baf+dafj; bgl = bgl+dbglj; bal = bal+dbalj;
        Rlf0 = Rlf0*LieExp(dqi); tlf0 = tlf0+dti; vlf0 = vlf0+dvi;
        bgf0 = bgf0+dbgfi; baf0 = baf0+dbafi; bgl0 = bgl0+dbgli; bal0 = bal0+dbali;

        [del_a_f,del_v_f,del_p_f] = GetDeltaIntgrateMeasurement(del_a_f,del_v_f,del_p_f,JacoInt,dbgfi,dbafi);
        [del_a_l,del_v_l,del_p_l] = GetDeltaIntgrateMeasurement(del_a_l,del_v_l,del_p_l,JacoIntL,dbgli,dbali);
    end
    
    Rlf0 = Rlf; tlf0 = tlf; vlf0 = vlf; bgf0 = bgf; baf0 = baf; bgl0 = bgl; bal0 = bal;
    Info = Jacobian'*infoAll*Jacobian;
    H11 = Info(1:indState,1:indState); 
    H10 = Info(1:indState,indState+1:2*indState); 
    H01 = Info(indState+1:2*indState,1:indState); 
    H00 = Info(indState+1:2*indState,indState+1:2*indState);
    infoPrior0 = H11-H10*inv(H00)*H01;
    x_est(i,:) = [rotm2quat(Rlf), tlf', vlf', bgf', baf', bgl', bal'];
    p_cov(i,:) = diag(inv(infoPrior0));
end
    
end

