function qhnew = RK4timestep(qh, dt)
    global filter
    k1 = RHS_eq(qh);
    k2 = RHS_eq(qh+k1*dt/2);
    k3 = RHS_eq(qh+k2*dt/2);
    k4 = RHS_eq(qh+k3*dt);
    qhnew = filter.*(qh + (k1 + 2*k2 + 2*k3 + k4)*dt/6);
end