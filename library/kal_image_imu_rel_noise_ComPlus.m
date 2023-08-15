function [x_est,p_Cov]=kal_image_imu_rel_noise_ComPlus(IMUF, IMUL, imuParameter, imuParameterL, nf, Tcb, camParameters, bPixel, obsPixelPH, rPixel, x0, p0) 
% IMUL=IMUL_;
m=size(obsPixelPH,1);
x_est = zeros(m,16);
p_Cov = zeros(m,15);
qe=x0(1:4)'; pe=x0(5:7)'; ve=x0(8:10)'; bge=x0(11:13)'; bae=x0(14:16)';
x_est(1,:) = x0(1:16);
a_f2l=quat2rotm(qe');

% Covariance
sigvw = imuParameter(1);
sigva = imuParameter(2);
siguw = imuParameter(3);
sigua = imuParameter(4);
sigvwL = imuParameterL(1);
sigvaL = imuParameterL(2); 

pCov=p0(1:15,1:15);
p_Cov(1,:)=diag(pCov)';
qCov=zeros(21,21);
qCov(1:3,1:3)=sigvw^2*eye(3);
qCov(4:6,4:6)=siguw^2*eye(3);
qCov(7:9,7:9)=sigva^2*eye(3);
qCov(10:12,10:12)=sigua^2*eye(3);
qCov(13:15,13:15)=sigvwL^2*eye(3);
qCov(16:18,16:18)=sigvaL^2*eye(3);
qCov(19:21,19:21)=sigvwL^2*eye(3);

tPre = obsPixelPH(1,1);
[imuF, iteIMUF, imuL, iteIMUL] = GetIMUAtFirstObservation(tPre, IMUF, IMUL);
mpImuFromLastFrame = [tPre imuF(2:end) imuL(2:end)];
nFace = size(bPixel,2);
nTagPerFace = 4;
for i=2:m

    tCur = obsPixelPH(i,1);

    j=iteIMUF;
    while j<size(IMUF,1) && IMUF(j,1)<tCur
        j=j+1;
    end
    imuF = IMUF(iteIMUF:j,:);
    iteIMUF = j; % >tCur  
    
    j=iteIMUL;
    while j<size(IMUL,1) && IMUL(j,1)<tCur
        j=j+1;
    end
    imuL = IMUL(iteIMUL:j,:);
    iteIMUL = j;   
 
	n = size(imuF,1);
    m = size(imuL,1);
    if n==0 || m==0
        error("Empty IMU measurements vector!!!");
    end
    
    mpImuFromLastFrame = [mpImuFromLastFrame; ImuSynchronization(imuF, imuL, tCur)];
    nImu = size(mpImuFromLastFrame,1);
    for j=1:nImu-1        
        dt = mpImuFromLastFrame(j+1,1)-mpImuFromLastFrame(j,1);
		wge = mpImuFromLastFrame(j,2:4)'-bge;
		acce = mpImuFromLastFrame(j,5:7)'-bae;
		wgel = mpImuFromLastFrame(j,8:10)';
		accel = mpImuFromLastFrame(j,11:13)';
        wgel2 = mpImuFromLastFrame(j+1,8:10)';

        % Partials
        a_fup=LieExp(dt*wge');
        a_lup=LieExp(dt*wgel');
        Jr=GetJacobian(dt*wge',2);
        JrL=GetJacobian(dt*wgel',2);
        JrL2=GetJacobian(dt*wgel2',2);

        f11=a_fup';
        f12=-Jr*dt;
        g11=-Jr*dt;
        g12=(a_lup'*a_f2l*a_fup)'*JrL*dt;

        f21=-dt^2/2*a_lup'*a_f2l*GetCrossMat(acce);
        f22=a_lup'*(eye(3)+GetCrossMat(wgel)*dt);
        f23=a_lup'*dt;
        f24=-dt^2/2*a_lup'*a_f2l;
        g21=-dt^2/2*a_lup'*a_f2l; 
%         g22=dt^2/2*a_lup';
%         g23=-GetCrossMat(a_lup'*((eye(3)+GetCrossMat(wgel)*dt)*pe+ve*dt+dt^2/2*(a_f2l*acce-accel)))*JrL*dt+...
%             a_lup'*GetCrossMat(pe)*dt;
        g23=dt^2/2*a_lup';
        g22=-GetCrossMat(a_lup'*((eye(3)+GetCrossMat(wgel)*dt)*pe+ve*dt+dt^2/2*(a_f2l*acce-accel)))*JrL*dt+...
            a_lup'*GetCrossMat(pe)*dt;
        
        f31=(dt/2*GetCrossMat(wgel2)-eye(3))*a_lup'*a_f2l*GetCrossMat(acce)*dt;
        f32=(GetCrossMat(wgel)-GetCrossMat(wgel2)-GetCrossMat(wgel2)*GetCrossMat(wgel)*dt)*a_lup';
        f33=(eye(3)-GetCrossMat(wgel2)*dt)*a_lup';
        f34=(dt/2*GetCrossMat(wgel2)-eye(3))*a_lup'*a_f2l*dt;
        g31=(dt/2*GetCrossMat(wgel2)-eye(3))*a_lup'*a_f2l*dt;
        g32=(eye(3)-GetCrossMat(wgel2)*dt)*(GetCrossMat(a_lup'*pe)-GetCrossMat(a_lup'*ve)*JrL2*dt)-...
            (GetCrossMat(wgel)-GetCrossMat(wgel2)-GetCrossMat(wgel2)*GetCrossMat(wgel)*dt)*GetCrossMat(a_lup'*pe)*JrL*dt-...
            (eye(3)+dt/2*GetCrossMat(wgel2))*GetCrossMat(a_lup'*(a_f2l*acce-accel))*JrL*dt*dt;
        g33=(eye(3)-dt/2*GetCrossMat(wgel2))*a_lup'*dt;
        g34=-GetCrossMat(a_lup'*ve)*dt-GetCrossMat((eye(3)+GetCrossMat(wgel)*dt)*a_lup'*pe)-...
            GetCrossMat(a_lup'*(a_f2l*acce-accel))*dt*dt/2;

        % Covariance Propagation
        f_mat=[f11 zeros(3,6) f12 zeros(3,3);
            f21 f22 f23 zeros(3,3) f24;
            f31 f32 f33 zeros(3,3) f34;
            zeros(3,9) eye(3) zeros(3,3);
            zeros(3,12) eye(3)];
        g_mat=[g11 zeros(3,9) g12 zeros(3,6); 
            zeros(3,6) g21 zeros(3,3) 1*g22 1*g23 zeros(3,3);
            zeros(3,6) 1*g31 zeros(3,3) 1*g32 1*g33 1*g34;
            zeros(3,3) dt*eye(3) zeros(3,15);
            zeros(3,9) dt*eye(3) zeros(3,9)];

        pCov=f_mat*pCov*f_mat'+g_mat*qCov*g_mat';

        % State Propagation
        a_f2l_=a_lup'*a_f2l*a_fup;
        pe_=a_lup'*((eye(3)+dt*GetCrossMat(wgel))*pe+ve*dt+0.5*dt^2*(a_f2l*acce-accel));
        ve_=(eye(3)-GetCrossMat(wgel2)*dt)*a_lup'*ve+(GetCrossMat(wgel)-GetCrossMat(wgel2)-GetCrossMat(wgel2)*GetCrossMat(wgel)*dt)*a_lup'*pe+...
            dt*(eye(3)-0.5*dt*GetCrossMat(wgel2))*a_lup'*(a_f2l*acce-accel);
        
        a_f2l = a_f2l_;
        pe = pe_;
        ve = ve_;
    end
    mpImuFromLastFrame = mpImuFromLastFrame(end,:);
    
    % Update
    r_est = []; h = []; obs=[];
    for j = 1:nFace
        if bPixel(i, j)==1
            rEstTemp = zeros(2*nTagPerFace,1);
            hTemp = zeros(2*nTagPerFace,15);
            for k = 1:nTagPerFace
                v3D = Tcb(1:3,:)*[a_f2l*nf(:,(j-1)*nTagPerFace+k)+pe; 1];
                rEstTemp((k-1)*2+1:k*2) = ProjectPinHole(v3D, camParameters);
                Jacobian = ProjectJacPinHole(v3D, camParameters);             
                hTemp((k-1)*2+1:k*2,1:3) = -Jacobian*Tcb(1:3,1:3)*a_f2l*GetCrossMat(nf(:,(j-1)*nTagPerFace+k));
                hTemp((k-1)*2+1:k*2,4:6) = Jacobian*Tcb(1:3,1:3);
            end
            r_est = [r_est; rEstTemp];
            obs = [obs; obsPixelPH(i,(j-1)*nTagPerFace*2+2:j*nTagPerFace*2+1)'];
            h = [h; hTemp];
        end
    end
    rcov = rPixel*eye(length(r_est));
    
    if ~isempty(r_est)
        gain=pCov*h'/(h*pCov*h'+rcov);
        xup=gain*(obs-r_est);    
    
        a_f2l=a_f2l*LieExp(xup(1:3)');  
        pe=pe+xup(4:6);
        ve=ve+xup(7:9);
        bge=bge+xup(10:12);
        bae=bae+xup(13:15);
        pCov=(eye(15)-gain*h)*pCov*(eye(15)-gain*h)'+gain*rcov*gain';
    else
        a_f2l=a_f2l;  
        pe=pe;
        ve=ve;
    end
    x_est(i,:) = [rotm2quat(a_f2l) pe' ve' bge' bae'];
    p_Cov(i,:)=diag(pCov)';
    
end

end