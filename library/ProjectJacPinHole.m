function [Jacobian] = ProjectJacPinHole(v3D, mvParameters)

    Jacobian = zeros(2,3);
    Jacobian(1, 1) = mvParameters(1) / v3D(3);
    Jacobian(1, 2) = 0;
    Jacobian(1, 3) = -mvParameters(1) * v3D(1) / (v3D(3) * v3D(3));
    Jacobian(2, 1) = 0;
    Jacobian(2, 2) = mvParameters(2) / v3D(3);
    Jacobian(2, 3) = -mvParameters(2) * v3D(2) / (v3D(3) * v3D(3));
end

