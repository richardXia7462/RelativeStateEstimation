function [errBiasGyro,jGyro,errBiasAcc,jAcc] = BiasErrJac(bg,bg0,ba,ba0)
    errBiasGyro = bg-bg0; 
    jacobianBiasVG = zeros(3,3);
    jacobianBiasVG(1:3,1:3) = eye(3);
    jacobianBiasVGk = zeros(3,3);
    jacobianBiasVGk(1:3,1:3) = -eye(3);
    jGyro = [zeros(3,9), jacobianBiasVG, zeros(3,12), jacobianBiasVGk, zeros(3,3)];

    errBiasAcc = ba-ba0; 
    jacobianBiasVA = zeros(3,3);
    jacobianBiasVA(1:3,1:3) = eye(3);
    jacobianBiasVAk = zeros(3,3);
    jacobianBiasVAk(1:3,1:3) = -eye(3);
    jAcc = [zeros(3,12), jacobianBiasVA, zeros(3,12), jacobianBiasVAk];
end