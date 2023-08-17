function [estOpt, diagP, residual] = GNSolver(iniEst, obsPixelPH, rPixel, bPixel, markF, Tcb, camParameters)

    nFace = size(bPixel,2);
    nTagPerFace = (size(obsPixelPH,2)-1)/nFace/2;
    num = size(obsPixelPH,1);
    estOpt = zeros(num,7);
    diagP = zeros(num,6);
    Rlf = quat2rotm(iniEst(1:4));
    tlf = iniEst(5:7)';
    threshold = 1e-5;
    residual = [];
    for i = 1:num
        obs=[];
        for j = 1:nFace
            if bPixel(i, j)==1
                obs = [obs; obsPixelPH(i,(j-1)*nTagPerFace*2+2:j*nTagPerFace*2+1)'];
            end
        end
        if isempty(obs)
            if i==1
                error("No feature at first frame!");
            elseif i==2
                estOpt(i,:)=estOpt(i-1,:);
                diagP(i,:)=diagP(i-1,:);   
            elseif size(estOpt,1)>2
                dtLast = obsPixelPH(i-1,1)-obsPixelPH(i-2,1);
                dtCur = obsPixelPH(i,1)-obsPixelPH(i-1,1);
                R1 = quat2rotm(estOpt(i-1,1:4)); R0 = quat2rotm(estOpt(i-2,1:4));
                t1 = estOpt(i-1,5:7); t0 = estOpt(i-2,5:7);
                Rlf = R1*LieExp(x_LieLog(R0\R1)/dtLast*dtCur);
                tlf = ((t1-t0)/dtLast*dtCur+t1)';
                estOpt(i,:) = [rotm2quat(Rlf) tlf'];
                diagP(i,:)=diagP(i-1,:); 
            end
            continue
        end
        rcov = rPixel*eye(length(obs));
        
        rep=1;
        ite=0;
        while rep==1 && ite<2
            v_est = []; h = []; 
            for j = 1:nFace
                if bPixel(i, j)==1
                    rEstTemp = zeros(2*nTagPerFace,1);
                    hTemp = zeros(2*nTagPerFace,6);
                    for k = 1:nTagPerFace
                        v3D = Tcb(1:3,:)*[Rlf*markF(:,(j-1)*nTagPerFace+k)+tlf; 1];
                        rEstTemp((k-1)*2+1:k*2) = ProjectPinHole(v3D, camParameters);
                        Jacobian = ProjectJacPinHole(v3D, camParameters);             
                        hTemp((k-1)*2+1:k*2,1:3) = -Jacobian*Tcb(1:3,1:3)*Rlf*GetCrossMat(markF(:,(j-1)*nTagPerFace+k));
                        hTemp((k-1)*2+1:k*2,4:6) = Jacobian*Tcb(1:3,1:3);
                    end
                    v_est = [v_est; rEstTemp];
                    h = [h; hTemp];
                end
            end

            K=inv(h'/rcov*h); 
            x=K*h'/rcov*(obs-v_est); % x=(J'/R*J)\J'/R*dR; 
%             norm(x)
%             norm(obs-r_est)
            if norm(x) > threshold
                axang=[x(1:3)'/norm(x(1:3)) norm(x(1:3))]; dt=x(4:6);
                Rlf=Rlf*axang2rotm(axang); 
                tlf=dt+tlf;
            else
                rep=0;
                estOpt(i,1:4)=rotm2quat(Rlf); 
                estOpt(i,5:7)=tlf';
                diagP(i,:)=diag(K);
                residual = [residual; norm(x)];
            end   
        end
    end
    
end

