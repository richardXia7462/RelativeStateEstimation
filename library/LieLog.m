function alg=LieLog(M)
	l=length(M);
    if l==4
        R=M(1:3,1:3);
        p=M(1:3,4);
        axang=rotm2axang(R);
        J=x_getJacobian(axang,1);
        r=J\p; phi=axang(1:3)*axang(4);
        alg=[r' phi];
    elseif l==3
        R=M(1:3,1:3);
        axang=rotm2axang(R);
        phi=axang(1:3)*axang(4);
        alg=phi;        
    else
        error('The size of the input is incorrect.')
    end
end