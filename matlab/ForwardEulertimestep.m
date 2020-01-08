function [qhnew] = ForwardEulertimestep(qh, dt)
    global filter
    RHS = RHS_eq(qh);    
    qhnew = filter.*(qh + RHS*dt);
end