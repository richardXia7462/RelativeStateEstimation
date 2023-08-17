function [del_a,del_v,del_p,Cov,Jacobian,CovG,CovA,correlation]=PreIntegrateMeasurement(Acc,Gyro,t,imuParameter,bg,ba)

    
    n=size(t,1);
    del_a = eye(3); del_v = zeros(3,1); del_p = zeros(3,1);
    JRg = zeros(3,3); JPa = zeros(3,3); JPg = zeros(3,3); JVa = zeros(3,3); JVg = zeros(3,3); 
    CovIMU = zeros(6);  
    CovIMU(1:3,1:3) = imuParameter(1)^2*eye(3); CovIMU(4:6,4:6) = imuParameter(2)^2*eye(3); 
    CovGyroWalk = imuParameter(3)^2*eye(3); 
    CovAccWalk = imuParameter(4)^2*eye(3);
    Cov = zeros(9); CovG = zeros(3); CovA = zeros(3); 
    del_a_J_E_0 = zeros(3,1);del_v_0 = zeros(3,1);del_p_0 = zeros(3,1);
    for i=1:n-1
        acc = Acc(i,:)'-ba; gyro = Gyro(i,:)'-bg;
        dt = t(i+1)-t(i);
        del_p = del_p+del_v*dt+0.5*del_a*acc*dt^2;
        del_v = del_v+del_a*acc*dt;
        del_a = del_a*LieExp(gyro'*dt);
        
        Jr = GetJacobian(gyro*dt,2);
        JPa = JPa+JVa*dt-0.5*JRg*dt^2;
        JPg = JPg+JVg*dt-0.5*JRg*dt^2*GetCrossMat(acc)*JRg;
        JVa = JVa-del_a*dt;
        JVg = JVg-del_a*dt*GetCrossMat(acc)*JRg;
        JRg = LieExp(gyro'*dt)'*JRg-Jr*dt;  
        
        A = [LieExp(gyro'*dt)', zeros(3), zeros(3);
            -del_a*GetCrossMat(acc)*dt, eye(3), zeros(3);
            -0.5*del_a*GetCrossMat(acc)*dt*dt, dt*eye(3), eye(3)];
        B = [Jr*dt, zeros(3);
            zeros(3), del_a*dt;
            zeros(3), 0.5*del_a*dt*dt]; 
        Cov = A*Cov*A'+B*CovIMU*B';
        CovG = CovG+CovGyroWalk; CovA = CovA+CovAccWalk;
        if i==1
            del_a_J_E_0 = del_a*Jr*CovIMU(1:3,1:3)*dt;
            del_v_0 = del_v;
            del_p_0 = del_p;
        end
    end
    dt = sum(diff(t))/(n-1);
    correlation = [del_a'*del_a_J_E_0; GetCrossMat(del_v_0-del_v)*del_a_J_E_0;...
        GetCrossMat(del_p_0-del_p)*del_a_J_E_0+(n-1)*GetCrossMat(del_v_0)*del_a_J_E_0*dt];
    Jacobian = [JRg,JVg,JVa,JPg,JPa];
end