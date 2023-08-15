%% load robot toolbox
% cd ./rvctools_1
% startup_rvc
% cd ..

%% configuration
cd library
clear, clc, close all
dt=0.001; 

% 0: disable 1: enable
bSynchroError = 0; 
bFullSmoother = 1; 
bNoisyCoordinate = 0;

%%
% parameters configuration
% Camera intrinsic parameters 
fx = 252.8263;
fy = 252.9636;
cx = 321.8137;
cy = 225.1054;
camParameters=[fx fy cx cy];
% Camera extrinsic parameters 
Tbc = [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1]'; % camera to body 
Tcb = inv(Tbc);
% Image size
width = 640;
height = 480;
range = [width height]'*ones(1,4);
sig3Pixel = 5;
        
tTag = 0.04;
kSam = tTag/dt;
% Detection threshold
thCosTheta = cos(pi*60/180);

imuRate = 200.0;
dtImu = 1/imuRate;
imuBatch = dtImu/dt; % 4ms an imu sample

Ng = 2.269e-03; Na = 8.182e-03;
Ngw = 1.536e-05; Naw = 6.154e-04;  
paraIMU = [Ng, Ngw, Na, Naw, 0, 0];
imuStdErrG = sqrt(imuRate)*Ng; imuStdErrA = sqrt(imuRate)*Na;
imuStdErrGW = 1/sqrt(imuRate)*Ngw; imuStdErrAW = 1/sqrt(imuRate)*Naw;
imuParameter=[imuStdErrG,imuStdErrA,imuStdErrGW,imuStdErrAW];

Ng = 1.528e-03; Na = 1.244e-03;
Ngw = 1.876e-05; Naw = 7.841e-04; 
paraIMUL = [Ng, Ngw, Na, Naw, 0, 0];
imuStdErrG = sqrt(imuRate)*Ng; imuStdErrA = sqrt(imuRate)*Na;
imuStdErrGW = 1/sqrt(imuRate)*Ngw; imuStdErrAW = 1/sqrt(imuRate)*Naw;
imuParameterL=[imuStdErrG,imuStdErrA,imuStdErrGW,imuStdErrAW];
sSynErr = bSynchroError*0.02; % synchronization error

% trajectory
radius = 0.2;
Height = 0.1;
kPer = 6;
acc = 9;

%%
centre = [0 0.5 -0.4];
kTraj = 6;
[AllTrajF, AllTrajL] = Trajectory(acc, kPer, Height, radius, dt, kTraj);

nS=size(AllTrajL,3);
% begin at t=0
tEnd=dt*(nS-1); 
tMov=(0:dt*1000:tEnd*1000)'/1000; 

quatF = tform2quat(AllTrajF); quatL = tform2quat(AllTrajL);
transF = tform2trvec(AllTrajF); transL = tform2trvec(AllTrajL);
velF=[diff(transF)/dt; zeros(1,3)]; velL=[diff(transL)/dt; zeros(1,3)];
accF=[diff(velF)/dt; zeros(1,3)]; accL=[diff(velL)/dt; zeros(1,3)];
velNormF = zeros(nS,1);
for i=1:nS 
    velNormF(i) = norm(velF(i,:));
end

figure, hold on, grid on
plot(tMov(1:nS/kPer/4),velNormF(1:nS/kPer/4))
xlabel('Time (sec)'), ylabel('Velocity (m/s)'), set(gca,'fontsize',28)

angRateF=zeros(nS,3); angRateL=zeros(nS,3); % angRateF: w^i_f|i
for i=1:nS-1   
    t1=AllTrajF(:,:,i); t2=AllTrajF(:,:,i+1); dt21=t2/t1;  
    angRateF(i,:)=LieLog(tform2rotm(dt21))/dt;
    t1=AllTrajL(:,:,i); t2=AllTrajL(:,:,i+1); dt21=t2/t1; 
    angRateL(i,:)=LieLog(tform2rotm(dt21))/dt;
end

% relative pose
quatF2L=zeros(nS,4); transFinL=zeros(nS,3); 
for i=1:nS
    Tf2l=AllTrajL(:,:,i)\AllTrajF(:,:,i);
    quatF2L(i,:)=tform2quat(Tf2l); transFinL(i,:)=tform2trvec(Tf2l);
end
velFinL=[diff(transFinL)/dt; zeros(1,3)];

figure
p1 = plot3(transF(:,1), transF(:,2), transF(:,3),'linewidth',1);
hold on, grid on, axis equal, set(gca,'fontsize',20)
p2 = plot3(transL(:,1), transL(:,2), transL(:,3),'linewidth',1);
Axis = 2*[0.1 0 0; 0 0.1 0; 0 0 0.1];    
Tf = AllTrajF(:,:,end);
temp = Tf(1:3,:)*[Axis; ones(1,3)];
plot3([Tf(1,4) temp(1,1)],[Tf(2,4) temp(2,1)],[Tf(3,4) temp(3,1)],'b','linewidth',1)
plot3([Tf(1,4) temp(1,2)],[Tf(2,4) temp(2,2)],[Tf(3,4) temp(3,2)],'r','linewidth',1)
plot3([Tf(1,4) temp(1,3)],[Tf(2,4) temp(2,3)],[Tf(3,4) temp(3,3)],'g','linewidth',1)
Tl = AllTrajL(:,:,end);
temp = Tl(1:3,:)*[Axis; ones(1,3)];
plot3([Tl(1,4) temp(1,1)],[Tl(2,4) temp(2,1)],[Tl(3,4) temp(3,1)],'b','linewidth',1)
plot3([Tl(1,4) temp(1,2)],[Tl(2,4) temp(2,2)],[Tl(3,4) temp(3,2)],'r','linewidth',1)
plot3([Tl(1,4) temp(1,3)],[Tl(2,4) temp(2,3)],[Tl(3,4) temp(3,3)],'g','linewidth',1)
legend([p1, p2], {'follower','leader'});
xlabel('x (m)'), ylabel('y (m)'), zlabel('z (m)')
% legend([p1, p2, p4], {'follower','leader','features'});
set(gca,'fontsize',24)

figure
p3 = plot(transFinL(:,1), transFinL(:,2),'g','linewidth',1);
% p3 = plot3(transFinL(:,1), transFinL(:,2), transFinL(:,3),'g','linewidth',1);
legend(p3, {'F|L'}); set(gca,'fontsize',34)
xlabel('x (m)'), ylabel('y (m)')
axis([-0.3 0.3 0.8 1.42])

imgSample = 1:kSam:nS;
nSamFea = length(imgSample);
quatF2LSam = quatF2L(imgSample,:);
transFinLSam = transFinL(imgSample,:);
velFinLSam = velFinL(imgSample,:);

imuSample = 1:imuBatch:nS;
nSimu = length(imuSample);
quatF2LIMU = quatF2L(imuSample,:);
transFinLIMU = transFinL(imuSample,:);
velFinLIMU = velFinL(imuSample,:);
tIMU = tMov(imuSample);
biasSam = tTag/dtImu;
biasSample = 1:biasSam:nSimu;

%% Marks Generation
tagsize = 0.14;
[nFace, nAncPerFac, markF, FNormal] = GetMarks(tagsize);
nF = size(markF,2); 
markFTruth = markF;
if bNoisyCoordinate
    markF = markF+0.005*randn(size(markF));
end

%% True visual feature measurement
Drop = [nSamFea+1, 8, 4];
kDrop = 1;
nDrop = Drop(kDrop);
[bPixel, obsPixelPHGT, plotTag] = GetVisObs(AllTrajL, AllTrajF, nFace, markFTruth, Tcb, thCosTheta, kSam, nSamFea, FNormal, camParameters, range, nDrop);
tSam = tMov(imgSample);
obsPixelPHGT = [tSam obsPixelPHGT];
% figure, hold on
% area(tSam/tTag, sum(bPixel,2),'DisplayName','bPixel')
% xlabel('Time (sec)'), 
% set(gca,'fontsize',20)

%% draw demo
bStep = 0;
% DrawTrajectoryObservation(transF, transL, AllTrajF, AllTrajL, imgSample, kSam, nFace, bPixel, plotTag, bStep);

%% add features noise
rPixel = (sig3Pixel/3)^2;
obsPixelPH = zeros(nSamFea,2*nF+1);
obsPixelPH(:,1) = obsPixelPHGT(:,1);
bFeaNoise = 1;
for i=1:size(bPixel,1)
    for j = 1:nFace
        if bPixel(i,j) == 1
            noise = bFeaNoise*normrnd(0,sig3Pixel/3,2,nAncPerFac);
            obsPixelPH(i,((j-1)*2*4+2:j*2*4+1)) = obsPixelPHGT(i,((j-1)*2*4+2:j*2*4+1))+noise(:)';
        end
    end
end

errAttNL=zeros(nSamFea,3,6); 
errTraNL=zeros(nSamFea,3,6); 
errVelNL=zeros(nSamFea,3,6); 
errBiasGNL=zeros(nSamFea,3,6); 
errBiasANL=zeros(nSamFea,3,6); 
estTrans = zeros(nSamFea,3,6); 

%% Nonlinear Optimization Using Extended Preintegration Full Biases
q0 = rotm2quat(quat2rotm(quatF2L(1,:))*axang2rotm([randn(1,3) 0.02]));
t0 = transFinL(1,:)+0.02*randn(1,3);
v0 = velFinL(1,:)+0.2*randn(1,3);
bg0 = 0.01*randn(1,3);
ba0 = 0.05*randn(1,3);
bgl0 = 0.01*randn(1,3);
bal0 = 0.05*randn(1,3);
iniEst = [q0 t0 v0 bg0 ba0 bgl0 bal0];
p0 = diag([0.02*ones(3,1); 0.02*ones(3,1); 0.2*ones(3,1); 0.01*ones(3,1); 0.05*ones(3,1); 0.01*ones(3,1); 0.05*ones(3,1);].^2);

%% ImgOnlyGN
'ImgOnlyGN'
[estOpt, diagP, residual] = GNSolver(iniEst, obsPixelPH, rPixel, bPixel, markF, Tcb, camParameters);
estTrans(:,:,5) = estOpt(:,5:7);
for i=1:nSamFea
    errAttNL(i,:,5)=LieLog(quat2rotm(estOpt(i,1:4))\quat2rotm(quatF2LSam(i,:)));
end
errTraNL(:,:,5)=transFinLSam-estOpt(:,5:7);

%% IMU simulation 
[accMeaL, gyroMeaL, ~, ~, bglTrue, balTrue] = GetIMU(accL(imuSample,:), angRateL(imuSample,:), quatL(imuSample,:), paraIMUL, dtImu, 1, 1);
IMUL = [tMov(imuSample) gyroMeaL accMeaL];
[accMeaF, gyroMeaF, accTrue, gyroTrue, bgfTrue, bafTrue] = GetIMU(accF(imuSample,:), angRateF(imuSample,:), quatF(imuSample,:), paraIMU, dtImu, 1, 1);
IMUF = [tMov(imuSample)+sSynErr gyroMeaF accMeaF];

biasGyroSam = bgfTrue(biasSample,:);
biasAccSam = bafTrue(biasSample,:);
biasGyroLSam = bglTrue(biasSample,:);
biasAccLSam = balTrue(biasSample,:);
GoundTruth = {quatF2LSam, transFinLSam, velFinLSam, biasGyroSam, biasAccSam, biasGyroLSam, biasAccLSam};

%% visual-inertial relative state estimation
% Full-Opt
'Full-Opt'
[x_est, p_cov] = NLOpt_GN_PH_ExtPreint_BothBias(IMUF, IMUL, iniEst, p0, obsPixelPH, rPixel, bPixel, markF, Tcb, camParameters, imuParameter, imuParameterL);
[errAttNL(:,:,4), errTraNL(:,:,4), errVelNL(:,:,4), errBiasGNL(:,:,4), errBiasANL(:,:,4)] = GetError(tSam, GoundTruth, x_est, p_cov, 0, 0, 0);
estTrans(:,:,4) = x_est(:,5:7);
PlotTwoBiases(obsPixelPH(:,1), [biasGyroSam biasGyroLSam], [x_est(:,11:13) x_est(:,17:19)], [p_cov(:,10:12) p_cov(:,16:18)])

% Minor-Opt
'Minor-Opt'
[x_est, p_cov] = extendedPreintegratedNonlinearOptimizationVis(IMUF, IMUL, imuParameter, imuParameterL, markF, Tcb, camParameters, bPixel, obsPixelPH, rPixel, iniEst, p0);
[errAttNL(:,:,1), errTraNL(:,:,1), errVelNL(:,:,1), errBiasGNL(:,:,1), errBiasANL(:,:,1)] = GetError(tSam, GoundTruth, x_est, p_cov, 1, 0, 0);
estTrans(:,:,1) = x_est(:,5:7);

% GBA-Opt
'GBA-Opt'
if bFullSmoother
    if bNoisyCoordinate
        [x_est, p_cov, mark_est] = extendedPreintegratedNonlinearOptimizationVisBA_mark(IMUF, IMUL, imuParameter, imuParameterL, markF, Tcb, camParameters, bPixel, obsPixelPH, rPixel, iniEst, p0);
    else
        [x_est, p_cov] = extendedPreintegratedNonlinearOptimizationVisBA(IMUF, IMUL, imuParameter, imuParameterL, markF, Tcb, camParameters, bPixel, obsPixelPH, rPixel, iniEst, p0);
    end
    [errAttNL(:,:,6), errTraNL(:,:,6), errVelNL(:,:,6), errBiasGNL(:,:,6), errBiasANL(:,:,6)] = GetError(tSam, GoundTruth, x_est, p_cov, 0, 0, 0);
    estTrans(:,:,6) = x_est(:,5:7);
end

% 2ndOrderEKF
'2ndOrderEKF'
[x_est, p_cov]=kal_image_imu_rel_noise_ComPlus(IMUF, IMUL, imuParameter, 1*imuParameterL, markF, Tcb, camParameters, bPixel, obsPixelPH, rPixel, iniEst, p0);
[errAttNL(:,:,2), errTraNL(:,:,2), errVelNL(:,:,2), errBiasGNL(:,:,2), errBiasANL(:,:,2)] = GetError(tSam, GoundTruth, x_est, p_cov, 0, 0, 0);
estTrans(:,:,2) = x_est(:,5:7);

% 3rdOrderEKF
'3rdOrderEKF'
[x_est, p_cov]=kal_image_imu_rel_noise_DiffPlus(IMUF, IMUL, imuParameter, imuParameterL, markF, Tcb, camParameters, bPixel, obsPixelPH, rPixel, iniEst, p0);
[errAttNL(:,:,3), errTraNL(:,:,3), errVelNL(:,:,3), errBiasGNL(:,:,3), errBiasANL(:,:,3)] = GetError(tSam, GoundTruth, x_est, p_cov, 0, 0, 0);
estTrans(:,:,3) = x_est(:,5:7);

%% Absolute trajectory error
rmse = zeros(1,6);
for i=1:nSamFea
    for j=1:6
        rmse(j) = rmse(j)+norm(errTraNL(i,:,j))^2;
    end
end
rmse = sqrt(rmse/nSamFea);
alg = [5,3,2,4,1,6];
Title={'ImaOnlyGN', '3rdOrderEKF', '2ndOrderEKF', 'FullOpt', 'MinorOpt', 'GBAOpt'};

figure
for i=1:6
    if i==6 && bFullSmoother==0
        continue
    end
    subplot(2,3,i), hold on, grid on
    p1=plot(estTrans(:,1,alg(i)), estTrans(:,2,alg(i)),'b','linewidth',1);
    p2=plot(transFinL(:,1), transFinL(:,2),'m','linewidth',3);
    xlabel('x (m)'), ylabel('y (m)')
    set(gca,'fontsize',20)
    title(Title{i})
    gtext(['RMSE: ',num2str(rmse(alg(i)))],'fontsize',15)
end
legend('Estimation','GroundTruth')

if bNoisyCoordinate && bFullSmoother
    errMark1 = GetErrorMark(markFTruth, markF);
    errMark2 = GetErrorMark(markFTruth, mark_est);
    figure, hold on, grid on
    plot([errMark1 errMark2],'LineWidth',3)
    set(gca,'fontsize',20)
    xlabel('feature ID')
    title('Distance to true coordinates')
    legend('raw','result')
end