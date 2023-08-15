function [bPixel, obsPixelPHGT, plotTag] = GetVisObs(AllTrajL, AllTrajF, nFace, markF, Tcb, thCosTheta, kSam, nSamFea, FNormal, camParameters, range, nDrop)

nS = size(AllTrajL,3);
nF = size(markF,2);
bPixel = zeros(nSamFea,nFace);
obsPixelPHGT = zeros(nSamFea,2*nF);
plotTag = zeros(nSamFea,3*nF);
for i=1:kSam:nS  
    
    if mod(round(i/kSam)+1, nDrop) == 0
        continue
    end
    
    bPixelTemp = zeros(1,nFace);
    obsPixelPHTemp = zeros(1,2*nF);
    plotTagTemp = zeros(1,3*nF);
    
    Tlf = AllTrajL(:,:,i)\AllTrajF(:,:,i);
    Sf = Tcb(1:3,:)*[Tlf(1:3,:)*[markF; ones(1,nF)]; ones(1,nF)];
    Sfi = AllTrajF(1:3,:,i)*[markF; ones(1,nF)];
    
    % up normal
    cosAngle = -(Tcb(1:3,1:3)*Tlf(1:3,1:3)*FNormal(:,1));
    cosAngle = cosAngle(3);
    if cosAngle > thCosTheta
    % Tag中心与光心连线和过Tag的平面法线夹角小于60度且相机坐标系下的z坐标为正
        obsPH = ProjectPinHole(Sf(:,1:4), camParameters);
        if sum(obsPH>range,'all')==0 && sum(obsPH<0,'all')==0
            bPixelTemp(1) = 1;
            obsPixelPHTemp(0*2*4+1:1*2*4) = obsPH(:);
            plotTagTemp(0*3*4+1:1*3*4) = reshape(Sfi(:,1:4),1,12);
        end
    end
    
    % front normal
    cosAngle = -(Tcb(1:3,1:3)*Tlf(1:3,1:3)*FNormal(:,2));
    cosAngle = cosAngle(3);
    if cosAngle > thCosTheta
    % Tag中心与光心连线和过Tag的平面法线夹角小于60度且相机坐标系下的z坐标为正
        obsPH = ProjectPinHole(Sf(:,5:8), camParameters);
        if sum(obsPH>range,'all')==0 && sum(obsPH<0,'all')==0
            bPixelTemp(2) = 1;
            obsPixelPHTemp(1*2*4+1:2*2*4) = obsPH(:);
            plotTagTemp(1*3*4+1:2*3*4) = reshape(Sfi(:,5:8),1,12);
        end
    end
    
    % back normal
    cosAngle = -(Tcb(1:3,1:3)*Tlf(1:3,1:3)*FNormal(:,3));
    cosAngle = cosAngle(3);
    if cosAngle > thCosTheta
    % Tag中心与光心连线和过Tag的平面法线夹角小于60度且相机坐标系下的z坐标为正
        obsPH = ProjectPinHole(Sf(:,9:12), camParameters);
        if sum(obsPH>range,'all')==0 && sum(obsPH<0,'all')==0
            bPixelTemp(3) = 1;
            obsPixelPHTemp(2*2*4+1:3*2*4) = obsPH(:);  
            plotTagTemp(2*3*4+1:3*3*4) = reshape(Sfi(:,9:12),1,12);
        end
    end
    
    % left normal
    cosAngle = -(Tcb(1:3,1:3)*Tlf(1:3,1:3)*FNormal(:,4));
    cosAngle = cosAngle(3);
    if cosAngle > thCosTheta
    % Tag中心与光心连线和过Tag的平面法线夹角小于60度且相机坐标系下的z坐标为正
        obsPH = ProjectPinHole(Sf(:,13:16), camParameters);
        if sum(obsPH>range,'all')==0 && sum(obsPH<0,'all')==0
            bPixelTemp(4) = 1;
            obsPixelPHTemp(3*2*4+1:4*2*4) = obsPH(:);
            plotTagTemp(3*3*4+1:4*3*4) = reshape(Sfi(:,13:16),1,12);
        end
    end
    
    % right normal   
    cosAngle = -(Tcb(1:3,1:3)*Tlf(1:3,1:3)*FNormal(:,5));
    cosAngle = cosAngle(3);
    if cosAngle > thCosTheta
    % Tag中心与光心连线和过Tag的平面法线夹角小于60度且相机坐标系下的z坐标为正
        obsPH = ProjectPinHole(Sf(:,17:20), camParameters);
        if sum(obsPH>range,'all')==0 && sum(obsPH<0,'all')==0
            bPixelTemp(5) = 1;
            obsPixelPHTemp(4*2*4+1:5*2*4) = obsPH(:);
            plotTagTemp(4*3*4+1:5*3*4) = reshape(Sfi(:,17:20),1,12);
        end
    end   
    
    % bottom normal   
    cosAngle = -(Tcb(1:3,1:3)*Tlf(1:3,1:3)*FNormal(:,6));
    cosAngle = cosAngle(3);
    if cosAngle > thCosTheta
    % Tag中心与光心连线和过Tag的平面法线夹角小于60度
        obsPH = ProjectPinHole(Sf(:,21:24), camParameters);
        if sum(obsPH>range,'all')==0 && sum(obsPH<0,'all')==0
            bPixelTemp(6) = 1;
            obsPixelPHTemp(5*2*4+1:6*2*4) = obsPH(:);
            plotTagTemp(5*3*4+1:6*3*4) = reshape(Sfi(:,21:24),1,12);
        end
    end 

    plotTag(round(i/kSam)+1, :) = plotTagTemp;
    bPixel(round(i/kSam)+1, :) = bPixelTemp;
    obsPixelPHGT(round(i/kSam)+1, :) = obsPixelPHTemp;
end
