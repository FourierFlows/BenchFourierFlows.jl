function [qhnew] = ForwardEulertimestep(qh, dt)
    RHS = RHS_eq(qh);    
    qhnew = qh + RHS*dt;
end