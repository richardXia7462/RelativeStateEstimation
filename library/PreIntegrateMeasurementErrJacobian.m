function jPreInt = PreIntegrateMeasurementErrJacobian(del_a_f,del_v_f,del_p_f,JacoInt,a_f2l0,a_f2l,del_a_l,gyrol0,gyrol1,T)

    Mr = del_a_l'*a_f2l0*del_a_f;
    dMr = LieLog(a_f2l'*Mr)';
    JRg = JacoInt(:,1:3); JVg = JacoInt(:,4:6); JVa = JacoInt(:,7:9);
    JPg = JacoInt(:,10:12); JPa = JacoInt(:,13:15);

    jacobianVP = zeros(9,6);
    invJl = inv(GetJacobian(dMr,1));
    jacobianVP(1:3,1:3) = -invJl;
    jacobianVP(4:6,4:6) = -del_a_l*GetCrossMat(gyrol1);
    jacobianVP(7:9,4:6) = -del_a_l;
    jacobianVV = zeros(9,3);
    jacobianVV(4:6,1:3) = -del_a_l;
    jacobianVPk = zeros(9,6);
    invJr = inv(GetJacobian(dMr,2));
    jacobianVPk(1:3,1:3) = invJr*del_a_f';
    jacobianVPk(4:6,1:3) = -a_f2l0*GetCrossMat(del_v_f);
    jacobianVPk(4:6,4:6) = GetCrossMat(gyrol0);
    jacobianVPk(7:9,1:3) = -a_f2l0*GetCrossMat(del_p_f);
    jacobianVPk(7:9,4:6) = eye(3)+GetCrossMat(gyrol0)*T;
    jacobianVVk = zeros(9,3);
    jacobianVVk(4:6,1:3) = eye(3);
    jacobianVVk(7:9,1:3) = eye(3)*T;
    jacobianVGk = zeros(9,3);
    jacobianVGk(1:3,1:3) = invJr*JRg;
    jacobianVGk(4:6,1:3) = a_f2l0*JVg;
    jacobianVGk(7:9,1:3) = a_f2l0*JPg;
    jacobianVAk = zeros(9,3);
    jacobianVAk(4:6,1:3) = a_f2l0*JVa;
    jacobianVAk(7:9,1:3) = a_f2l0*JPa;           
%     jPreInt = [jacobianVP, jacobianVV, jacobianVPk, jacobianVVk, jacobianVGk, jacobianVAk];   
    jPreInt = [jacobianVP, jacobianVV, zeros(9,6), jacobianVPk, jacobianVVk, jacobianVGk, jacobianVAk];    
end