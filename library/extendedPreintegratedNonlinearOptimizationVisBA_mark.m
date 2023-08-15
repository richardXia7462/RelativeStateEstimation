function [x_est, p_cov, mark_est] = extendedPreintegratedNonlinearOptimizationVisBA_mark(IMUF, IMUL, imuParameter, imuParameterL, markFtoEst, Tcb, camParameters, bPixel, obsPixelPH, rPixel, x0, p0)
% x0 = iniEst;

nS = size(obsPixelPH,1);
% nS = 20;
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

% global BA
EPI = cell(nS-1,12);
RW = cell(nS-1,2);
FO = cell(nS); 
STATE = cell(nS,5);
STATE(1,:) = {Rlf0, tlf0, vlf0, bg0, ba0};
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
    EPI(i-1,:) = {del_a_f, del_v_f, del_p_f, CovF, JacoInt, del_a_l, del_v_l, del_p_l, CovL, correl, gyrol0,gyrol1};
    RW(i-1,:) = {CovG, CovA};
    
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
    FO(i) = {obs};
    
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
                    v3D = Tcb(1:3,:)*[Rlf*markFtoEst(:,(j-1)*nTagPerFace+k)+tlf; 1];
                    rEstTemp((k-1)*2+1:k*2) = Project_PinHole(v3D, camParameters);
                    Jac = projectJac_PinHole(v3D, camParameters);             
                    hTemp((k-1)*2+1:k*2,1:3) = -Jac*Tcb(1:3,1:3)*Rlf*x_crossMat(markFtoEst(:,(j-1)*nTagPerFace+k));
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
%         Norm = norm(errAll+Jacobian*deltaX)
        
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
    
    STATE(i,:) = {Rlf0, tlf0, vlf0, bg0, ba0};
end

% First observation
obs=[];
for j = 1:nFace
    if bPixel(1, j)==1
        obs = [obs; obsPixelPH(1,(j-1)*nTagPerFace*2+2:j*nTagPerFace*2+1)'];
    end
end
FO{1} = obs;

mark_est = markFtoEst;
szMark = nTagPerFace*nFace*3;
for ite=1:3
    % First observation
    i = 1;
    Rlf = STATE{i,1}; tlf = STATE{i,2}; 
    v_est = []; h = []; hMark = [];
    for j = 1:nFace
        if bPixel(i, j)==1
            rEstTemp = zeros(2*nTagPerFace, 1);
            hTemp = zeros(2*nTagPerFace, 15);
            hMarkTemp = zeros(2*nTagPerFace, szMark);
            iFirstMark = (j-1)*nTagPerFace*3;
            for k = 1:nTagPerFace
                v3D = Tcb(1:3,:)*[Rlf*mark_est(:,(j-1)*nTagPerFace+k)+tlf; 1];
                rEstTemp((k-1)*2+1:k*2) = ProjectPinHole(v3D, camParameters);
                Jac = ProjectJacPinHole(v3D, camParameters);             
                hTemp((k-1)*2+1:k*2,1:3) = -Jac*Tcb(1:3,1:3)*Rlf*x_crossMat(mark_est(:,(j-1)*nTagPerFace+k));
                hTemp((k-1)*2+1:k*2,4:6) = Jac*Tcb(1:3,1:3);
                % mark jacobian
                hMarkTemp((k-1)*2+1:k*2,iFirstMark+((k-1)*3+1:k*3))=Jac*Tcb(1:3,1:3)*Rlf;
            end
            v_est = [v_est; rEstTemp];
            h = [h; hTemp];
            hMark = [hMark; hMarkTemp];
        end
    end
    if ~isempty(h)
        h = [h(:,16:end) h(:,1:15)];
    else
        error('No first visual measurement')
    end
    
    errFeature = FO{1}-v_est;
    jacFeature = -h;
    jacMark = -hMark;
    infoVis = 1/rPixel*eye(length(FO{1}));

    Info = infoVis;
    errAll = errFeature;

    Jacobian = [jacMark jacFeature];
    
    for i=1:nS-1
        del_a_f = EPI{i,1}; del_v_f = EPI{i,2}; del_p_f = EPI{i,3};
        del_a_l = EPI{i,6}; del_v_l = EPI{i,7}; del_p_l = EPI{i,8};
        gyrol0 = EPI{i,11}; gyrol1 = EPI{i,12};
        Rlf0 = STATE{i,1}; tlf0 = STATE{i,2}; vlf0 = STATE{i,3}; bg0 = STATE{i,4}; ba0 = STATE{i,5}; 
        Rlf = STATE{i+1,1}; tlf = STATE{i+1,2}; vlf = STATE{i+1,3}; bg = STATE{i+1,4}; ba = STATE{i+1,5}; 
        tij = obsPixelPH(i+1,1)-obsPixelPH(i,1);
        Mr = del_a_l'*Rlf0*del_a_f;
        dMr = x_LieLog(Rlf'*Mr)';
        dMv = Rlf0*del_v_f-del_v_l+vlf0+x_crossMat(gyrol0)*tlf0-del_a_l*(vlf+x_crossMat(gyrol1)*tlf);
        dMp = Rlf0*del_p_f-del_p_l-del_a_l*tlf+tlf0+(vlf0+x_crossMat(gyrol0)*tlf0)*tij;
        errPreIntegration = [dMr; dMv; dMp];

        CovF = EPI{i,4}; CovL = EPI{i,9}; correl = EPI{i,10};
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

        CovG = RW{i,1}; CovA = RW{i,2};
        infoBiasGyro = inv(CovG); infoBiasAcc = inv(CovA);

        JacoInt = EPI{i,5};
        jPreInt = preIntegrateMeasurementErrJacobian(del_a_f,del_v_f,del_p_f,JacoInt,Rlf0,Rlf,del_a_l,gyrol0,gyrol1,tij);   
        jPreInt = [jPreInt(:,16:end) jPreInt(:,1:15)];
        
        v_est = []; h = []; hMark = [];
        for j = 1:nFace
            if bPixel(i+1, j)==1
                rEstTemp = zeros(2*nTagPerFace,1);
                hTemp = zeros(2*nTagPerFace,30);
                hMarkTemp = zeros(2*nTagPerFace, szMark);
                iFirstMark = (j-1)*nTagPerFace*3;
                for k = 1:nTagPerFace
                    v3D = Tcb(1:3,:)*[Rlf*mark_est(:,(j-1)*nTagPerFace+k)+tlf; 1];
                    rEstTemp((k-1)*2+1:k*2) = Project_PinHole(v3D, camParameters);
                    Jac = projectJac_PinHole(v3D, camParameters);             
                    hTemp((k-1)*2+1:k*2,1:3) = -Jac*Tcb(1:3,1:3)*Rlf*x_crossMat(mark_est(:,(j-1)*nTagPerFace+k));
                    hTemp((k-1)*2+1:k*2,4:6) = Jac*Tcb(1:3,1:3);
                    % mark jacobian
                    hMarkTemp((k-1)*2+1:k*2,iFirstMark+((k-1)*3+1:k*3))=Jac*Tcb(1:3,1:3)*Rlf;
                end
                v_est = [v_est; rEstTemp];
                h = [h; hTemp];
                hMark = [hMark; hMarkTemp];
            end
        end
        if ~isempty(h)
            h = [h(:,16:end) h(:,1:15)];
        end
        obs = FO{i+1};
        errFeature = obs-v_est;
        jacFeature = -h;
        jacMark = -hMark;
        infoVis = 1/rPixel*eye(length(obs));

        [errBiasGyro,jGyro,errBiasAcc,jAcc] = biasErrJac(bg,bg0,ba,ba0);
        jGyro = [jGyro(:,16:end) jGyro(:,1:15)];
        jAcc = [jAcc(:,16:end) jAcc(:,1:15)];
        
        tempErr = [errPreIntegration; errBiasGyro; errBiasAcc; errFeature];
        tempInfo = zeros(15+length(obs), 15+length(obs));
        tempInfo(1:9,1:9) = infoPreIntegration;
        tempInfo(10:12,10:12) = infoBiasGyro;
        tempInfo(13:15,13:15) = infoBiasAcc;
        tempInfo(16:end,16:end) = infoVis;

        Info = [Info zeros(size(errAll,1),size(tempErr,1)); zeros(size(errAll,1),size(tempErr,1))' tempInfo];
        errAll = [errAll; tempErr];
        tempJac = [jPreInt; jGyro; jAcc; jacFeature];
%         Jacobian = [Jacobian zeros(size(Jacobian,1),15); zeros(size(tempErr,1),size(Jacobian,2)-15) tempJac];
        jacMark = [zeros(15, szMark); jacMark];
        Jacobian = [Jacobian zeros(size(Jacobian,1),15); jacMark zeros(size(tempErr,1), size(Jacobian,2)-15-szMark) tempJac];
    end
    
    del = sum(bPixel)==0;
    for iDel=1:nFace
        if del(iDel)==1
            Jacobian(:,(iDel-1)*nTagPerFace*3:iDel*nTagPerFace*3) = [];
        end
    end
%     Norm = norm(errAll)
    errAll'*Info*errAll
    tic
    deltaX = -(Jacobian'*Info*Jacobian)\Jacobian'*Info*errAll; %norm(deltaX) 
    toc
    (errAll+Jacobian*deltaX)'*Info*(errAll+Jacobian*deltaX)
%     deltaX = -(Jacobian'*Jacobian)\Jacobian'*errAll; %norm(deltaX)
%     Norm = norm(errAll+Jacobian*deltaX)
    nDel = 0;
    for i=1:nFace
        if del(iDel)==1
            nDel = nDel+1;
        end
        iFirst = (i-1)*nTagPerFace*3;
        dMarki = reshape(deltaX(iFirst+(1:nTagPerFace*3)),3,nTagPerFace);
        mark_est(:,(i-1)*nTagPerFace+1:i*nTagPerFace) = mark_est(:,(i-1)*nTagPerFace+1:i*nTagPerFace)+dMarki;
    end
    for i=1:nS
        index = 15*(i-1)+szMark;
        dq = deltaX(index+(1:3)); dt = deltaX(index+(4:6)); dv = deltaX(index+(7:9)); dbg = deltaX(index+(10:12)); dba = deltaX(index+(13:15)); 
        Rlf = STATE{i,1}; tlf = STATE{i,2}; vlf = STATE{i,3}; bg = STATE{i,4}; ba = STATE{i,5}; 
        Rlf = Rlf*x_LieExp(dq); tlf = tlf+dt; vlf = vlf+dv; bg = bg+dbg; ba = ba+dba;
        STATE(i,:) = {Rlf, tlf, vlf, bg, ba};
        if i<nS
            JacoInt = EPI{i,5};
            del_a_f = EPI{i,1}; del_v_f = EPI{i,2}; del_p_f = EPI{i,3};
            [del_a_f,del_v_f,del_p_f] = getDeltaIntgrateMeasurement(del_a_f,del_v_f,del_p_f,JacoInt,dbg,dba);
            EPI{i,1} = del_a_f; EPI{i,2} = del_v_f; EPI{i,3} = del_p_f;
        end
    end
end

Cov = inv(Jacobian'*Info*Jacobian);
x_est = zeros(nS,16);
p_cov = zeros(nS,15);
for i=1:nS
    Rlf = STATE{i,1}; tlf = STATE{i,2}; vlf = STATE{i,3}; bg = STATE{i,4}; ba = STATE{i,5};
    x_est(i,:) = [rotm2quat(Rlf), tlf', vlf', bg', ba'];
    index = 15*(i-1)+szMark;
    p_cov(i,:) = diag(Cov(index+(1:15),index+(1:15)));
end

end

