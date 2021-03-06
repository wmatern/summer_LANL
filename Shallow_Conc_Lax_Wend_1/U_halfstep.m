function U = U_halfstep(deltaT, U_i, U_n, F_i, F_n, r_i, r_n, A_i,...
    A_n, V_i, V_n)
U = ((r_i*U_n+r_n*U_i)/(r_i+r_n))-0.5d0*deltaT*...
    ((F_n*A_n*min(1.0d0, A_i/A_n)-F_i*A_i*min(1.0d0, A_n/A_i))/...
    (V_n*min(0.5d0, V_i/V_n)+V_i*min(0.5d0, V_n/V_i)));