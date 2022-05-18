function dydt = ajm_1d(t, y, beta, Text, B, taum, tauv)

    dydt = zeros(3,1);
    Tst = 0.3;
    k0 = 2/Tst;
    u  = y(1); 
    m  = y(2);
    l0 = y(3);
    T = u + beta*(m - 0.5) + B*(u + l0 - 1);
    dydt(1) = -u/tauv - T + Text;
    dydt(2) = ( 1 - m*( 1 + exp(-k0*(T-Tst)) ) )/taum;
    dydt(3) = u/tauv;

end