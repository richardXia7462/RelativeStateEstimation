function T=LieExp(g)
    l=length(g);
    if l==6
        p=g(1:3);
        phi=g(4:6);
        if norm(phi)==0
            R=eye(3);
        else
            R=axang2rotm([phi/norm(phi) norm(phi)]);
        end
        r=x_getJacobian(phi)*p';
        T=[R r; zeros(1,3) 1];
    elseif l==3
        phi=g;
        if norm(phi)==0
            R=eye(3);
        else
            if(size(phi,1)==3)
                phi=phi';
            end
            R=axang2rotm([phi/norm(phi) norm(phi)]);
        end
        T=R;
    else
        error('The length of the argument is incorrect.')
    end

end