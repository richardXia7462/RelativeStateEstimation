function J=GetJacobian(liAlg,lr)
% lr: 1 left, 2 right
    size_ = size(liAlg);
    if(size_(1)>1 && size_(2)>1)
        error("Jacobian input must be vector ");
    end
    l=length(liAlg);
    if l==3 
        phi=norm(liAlg);
        a=liAlg/phi;
    elseif l==4
        a=liAlg(1:3)/norm(liAlg(1:3));
        phi=liAlg(4);
    end
    if norm(phi)==0
        J=eye(3);
    elseif lr==1
        J=sin(phi)/phi*eye(3)+(1-sin(phi)/phi)*(a*a')+(1-cos(phi))/phi*GetCrossMat(a);
    elseif lr==2
        J=sin(phi)/phi*eye(3)+(1-sin(phi)/phi)*(a*a')-(1-cos(phi))/phi*GetCrossMat(a);
    end
end