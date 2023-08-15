function [accm, wgm, acc, wg, gBias, aBias] = GetIMU(acci, wgi, quatB2I, paraIMU, dt, bBias, bNoise)

    m = size(quatB2I,1);
    gI = [0, 0, -9.81]';
    acc = zeros(m,3);
    wg = zeros(m,3);
    for i=1:m
        Ai2b = quat2rotm(quatB2I(i,:))';
        acc(i,:) = (Ai2b*(acci(i,:)'-gI))';
        wg(i,:) = (Ai2b*wgi(i,:)')';
    end
    sigv=paraIMU(1);sigu=paraIMU(2);
    sigva=paraIMU(3);sigua=paraIMU(4);
    ksg1=paraIMU(5);ksg2=paraIMU(5);ksg3=paraIMU(5);
    ksa1=paraIMU(6);ksa2=paraIMU(6);ksa3=paraIMU(6);
    
%     bNoise = 1;
    % Gyros
    num_g=dt*[1 1];den_g=2*[1 -1];
    [phi_g,gam_g,c_g,d_g]=tf2ss(num_g,den_g);
    gbias1=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),1*pi/180/3600);
    gbias2=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),1*pi/180/3600);
    gbias3=dlsim(phi_g,gam_g,c_g,d_g,sigu/sqrt(dt)*randn(m,1),1*pi/180/3600);

    gwhite1=sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1);
    gwhite2=sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1);
    gwhite3=sqrt(sigv^2/dt+1/12*sigu^2*dt)*randn(m,1);

    wgm1=(1+ksg1)*wg(:,1)+bBias*gbias1+bNoise*gwhite1;
    wgm2=(1+ksg2)*wg(:,2)+bBias*gbias2+bNoise*gwhite2;
    wgm3=(1+ksg3)*wg(:,3)+bBias*gbias3+bNoise*gwhite3;
    wgm=[wgm1 wgm2 wgm3];
    gBias = [gbias1 gbias2 gbias3];
    
    % Accelerometers
    num_a=dt*[1 1];den_a=2*[1 -1];
    [phi_a,gam_a,c_a,d_a]=tf2ss(num_a,den_a);
    abias1=dlsim(phi_a,gam_a,c_a,d_a,sigua/sqrt(dt)*randn(m,1),0.003);
    abias3=dlsim(phi_a,gam_a,c_a,d_a,sigua/sqrt(dt)*randn(m,1),0.003);
    abias2=dlsim(phi_a,gam_a,c_a,d_a,sigua/sqrt(dt)*randn(m,1),0.003);

    awhite1=sqrt(sigva^2/dt+1/12*sigua^2*dt)*randn(m,1);
    awhite2=sqrt(sigva^2/dt+1/12*sigua^2*dt)*randn(m,1);
    awhite3=sqrt(sigva^2/dt+1/12*sigua^2*dt)*randn(m,1);

    accm1=(1+ksa1)*acc(:,1)+bBias*abias1+bNoise*awhite1;
    accm2=(1+ksa2)*acc(:,2)+bBias*abias2+bNoise*awhite2;
    accm3=(1+ksa3)*acc(:,3)+bBias*abias3+bNoise*awhite3;
    accm=[accm1 accm2 accm3];
    aBias = [abias1 abias2 abias3];
end

