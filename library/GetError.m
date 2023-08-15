function [errAtt, errTra, errVel, errBiasG, errBiasA] = GetError(t, goundtruth, estimation, covariance, bPVQ, bBiasF, bBiasL)

    nSamFea = length(t);
    quatF2LSam = goundtruth{1};
    transFinLSam = goundtruth{2};
    velFinLSam = goundtruth{3};
    biasGyroSam = goundtruth{4};
    biasAccSam = goundtruth{5};

    errAtt = zeros(nSamFea, 3);
    for i=1:nSamFea
        errAtt(i,:)=LieLog(quat2rotm(estimation(i,1:4))\quat2rotm(quatF2LSam(i,:)));
    end
    errTra=transFinLSam-estimation(:,5:7);
    errVel=velFinLSam-estimation(:,8:10);
    errBiasG=biasGyroSam-estimation(:,11:13);
    errBiasA=biasAccSam-estimation(:,14:16);      
    if bBiasL
        biasGyroLSam = goundtruth{6};
        biasAccLSam = goundtruth{7};
        error = {errAtt, errTra, errVel, errBiasG, errBiasA, biasGyroLSam-estimation(:,17:19), biasAccLSam-estimation(:,20:22)};
        cov = {covariance(:,1:3), covariance(:,4:6), covariance(:,7:9), covariance(:,10:12), covariance(:,13:15), covariance(:,16:18), covariance(:,19:21)};
        DrawEstimation3Sigma(t, error, cov, bPVQ, bBiasF, bBiasL);
    else
        error = {errAtt, errTra, errVel, errBiasG, errBiasA};
        cov = {covariance(:,1:3), covariance(:,4:6), covariance(:,7:9), covariance(:,10:12), covariance(:,13:15)};
        DrawEstimation3Sigma(t, error, cov, bPVQ, bBiasF, bBiasL);
    end
