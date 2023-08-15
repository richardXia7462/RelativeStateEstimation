function [del_a,del_v,del_p] = GetDeltaIntgrateMeasurement(del_a,del_v,del_p,JacoInt,dbg,dba)
    JRg = JacoInt(:,1:3); JVg = JacoInt(:,4:6); JVa = JacoInt(:,7:9);
    JPg = JacoInt(:,10:12); JPa = JacoInt(:,13:15);
    del_a = del_a*LieExp(JRg*dbg);
    del_v = del_v+JVg*dbg+JVa*dba;
    del_p = del_p+JPg*dbg+JPa*dba;
end