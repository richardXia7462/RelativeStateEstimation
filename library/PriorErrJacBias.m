function [errPrior,jPrior] = PriorErrJacBias(priorState,estState)

    a_f2l0_mea = quat2rotm(priorState(1:4)'); p_f2l0_mea = priorState(5:7);
    v_f2l0_mea = priorState(8:10); bg0_mea = priorState(11:13); ba0_mea = priorState(14:16);
    bgl0_mea = priorState(17:19); bal0_mea = priorState(20:22);
    a_f2l0 = quat2rotm(estState(1:4)'); p_f2l0 = estState(5:7);
    v_f2l0 = estState(8:10); bg0 = estState(11:13); ba0 = estState(14:16);
    bgl0 = estState(17:19); bal0 = estState(20:22);
    
    erPrior = LieLog(a_f2l0_mea'*a_f2l0)';
    errPrior = [erPrior; p_f2l0-p_f2l0_mea; v_f2l0-v_f2l0_mea; 
        bg0-bg0_mea; ba0-ba0_mea; bgl0-bgl0_mea; bal0-bal0_mea];
    jacobianPriorVPk = zeros(21,6);
    invJrPrior = inv(GetJacobian(erPrior,2));
    jacobianPriorVPk(1:3,1:3) = invJrPrior;
    jacobianPriorVPk(4:6,4:6) = eye(3);
    jacobianPriorVVk = zeros(21,3);
    jacobianPriorVVk(7:9,1:3) = eye(3);
    jacobianPriorVGk = zeros(21,3);
    jacobianPriorVGk(10:12,1:3) = eye(3);
    jacobianPriorVAk = zeros(21,3);
    jacobianPriorVAk(13:15,1:3) = eye(3);
    jacobianPriorVGLk = zeros(21,3);
    jacobianPriorVGLk(16:18,1:3) = eye(3);
    jacobianPriorVALk = zeros(21,3);
    jacobianPriorVALk(19:21,1:3) = eye(3);
%     jPrior = [zeros(15,9), jacobianPriorVPk, jacobianPriorVVk, jacobianPriorVGk, jacobianPriorVAk];
    jPrior = [zeros(21,21), jacobianPriorVPk, jacobianPriorVVk, jacobianPriorVGk, jacobianPriorVAk, jacobianPriorVGLk, jacobianPriorVALk];


end