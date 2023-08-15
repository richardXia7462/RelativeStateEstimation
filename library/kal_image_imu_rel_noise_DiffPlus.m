function [x_est,p_cov]=kal_image_imu_rel_noise_DiffPlus(IMUF, IMUL, imuParameter, imuParameterL, nf, Tcb, camParameters, bPixel, obsPixelPH, rPixel, x0, p0)                 

NS=size(obsPixelPH,1);
x_est = zeros(NS,16);
p_cov = zeros(NS,15);

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
p_cov(1,:)=diag(pCov)';
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
for i=2:NS
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
        if dt==0
            continue
        end
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

        f11 = a_fup';
        f14 = -Jr*dt;
        g11 = -Jr*dt;
        g15 = (a_lup'*a_f2l*a_fup)'*JrL*dt;

        da1 = -a_f2l*GetCrossMat(acce);
        da2 = 1/dt*GetCrossMat(wgel-wgel2)-GetCrossMat(wgel)*GetCrossMat(wgel);
        da3 = -2*GetCrossMat(wgel);
        da5 = -a_f2l;
        da3_ = -a_f2l;
        da5_ = 1/dt*GetCrossMat(pe)-(GetCrossMat(wgel)*GetCrossMat(pe)+GetCrossMat(GetCrossMat(wgel)*pe))-2*GetCrossMat(ve);
        da6_ = eye(3);
        da7_ = -1/dt*GetCrossMat(pe);

        f21 = 0.5*dt*dt*da1;
        f22 = eye(3)+0.5*dt*dt*da2;
        f23 = dt*eye(3)+0.5*dt*dt*da3;
        f25 = 0.5*dt*dt*da5;
        g23 = 0.5*dt*dt*da3_;
        g25 = 0.5*dt*dt*da5_;
        g26 = 0.5*dt*dt*da6_;
        g27 = 0.5*dt*dt*da7_;

        f31 = dt*da1;
        f32 = dt*da2;
        f33 = eye(3)+dt*da3;
        f35 = dt*da5;
        g33 = dt*da3_;
        g35 = dt*da5_;
        g36 = dt*da6_;
        g37 = dt*da7_;

        % Covariance Propagation
        f_mat=[f11 zeros(3,6) f14 zeros(3,3);
               f21 f22 f23 zeros(3,3) f25;
               f31 f32 f33 zeros(3,3) f35;
               zeros(3,9) eye(3) zeros(3,3);
               zeros(3,12) eye(3)];
        g_mat=[g11 zeros(3,9) g15 zeros(3,6); 
               zeros(3,6) g23 zeros(3,3) g25 g26 g27;
               zeros(3,6) g33 zeros(3,3) g35 g36 g37;
               zeros(3,3) dt*eye(3) zeros(3,15);
               zeros(3,9) dt*eye(3) zeros(3,9)];
        
        pCov=f_mat*pCov*f_mat'+g_mat*qCov*g_mat';

        % State Propagation
        a_f2l_=a_lup'*a_f2l*a_fup;
        ae = -1/dt*GetCrossMat(wgel2-wgel)*pe-GetCrossMat(wgel)*GetCrossMat(wgel)*pe-2*GetCrossMat(wgel)*ve+a_f2l*acce-accel;
        ve_ = ve+ae*dt;
        pe_ = pe+ve*dt+0.5*ae*dt*dt;
        
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
                rEstTemp((k-1)*2+1:k*2) = Project_PinHole(v3D, camParameters);
                Jacobian = projectJac_PinHole(v3D, camParameters);             
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
    
        a_f2l=a_f2l*x_LieExp(xup(1:3)');  
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
    p_cov(i,:)=diag(pCov)';

end

end