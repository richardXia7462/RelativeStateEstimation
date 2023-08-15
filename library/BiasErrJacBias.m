function [errBiasGyro,jGyro,errBiasAcc,jAcc] = biasErrJacBias(bg, bg0, ba, ba0, tij)
    errBiasGyro = bg-bg0;
    errBiasGyro = errBiasGyro/tij;
    jacobianBiasVG = zeros(3,3);
    jacobianBiasVG(1:3,1:3) = 1/tij*eye(3);
    jacobianBiasVGk = zeros(3,3);
    jacobianBiasVGk(1:3,1:3) = -1/tij*eye(3);
    jGyro = {jacobianBiasVG, jacobianBiasVGk};
  
    errBiasAcc = ba-ba0;
    errBiasAcc = errBiasAcc/tij;
    jacobianBiasVA = zeros(3,3);
    jacobianBiasVA(1:3,1:3) = 1/tij*eye(3);
    jacobianBiasVAk = zeros(3,3);
    jacobianBiasVAk(1:3,1:3) = -1/tij*eye(3);
    jAcc = {jacobianBiasVA, jacobianBiasVAk};
end