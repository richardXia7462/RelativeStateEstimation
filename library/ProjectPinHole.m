function [p_pix] = ProjectPinHole(p_cam, mvParameters)

    n = size(p_cam,2);
    p_pix = zeros(2,n);
    for i=1:n
        p3D = p_cam(:,i);
        pix = [mvParameters(1) * p3D(1) / p3D(3) + mvParameters(3);
                mvParameters(2) * p3D(2) / p3D(3) + mvParameters(4)];
        p_pix(:,i) = pix;
    end
end

