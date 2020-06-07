function [qhnew, RHSp, RHSpp] = AB3timestep(qh, RHSp, RHSpp, dt)
    RHS = RHS_eq(qh);   
    qhnew = qh + (23*RHS - 16*RHSp + 5*RHSpp)*dt/12;
    RHSpp = RHSp;
    RHSp  = RHS;
end