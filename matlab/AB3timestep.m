function [qhnew, RHSp, RHSpp] = AB3timestep(qh, RHSp, RHSpp, dt)
    global filter
    RHS = RHS_eq(qh);   
    qhnew = filter.*(qh + (23*RHS - 16*RHSp + 5*RHSpp)*dt/12);
    RHSpp = RHSp;
    RHSp  = RHS;
end