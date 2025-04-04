a = 0.85;
a1 = a; a2 = a; a3 = a;
tc = 0.85; tw = 0.85;
tc_low = tc; tc_upp = tc;
phi = pi:0.01*pi:5*pi;
T = zeros(1, length(phi));
counter1 = 1;

for jj=phi
    ph1 = jj;
    ph2 = jj;
    ph3 = jj;

    z_upp = (tc_upp-a3*exp(-(1i)*ph3))/(1-tc_upp*a3*exp(-(1i)*ph3));
    z_low = (tc_low-a2*z_upp*exp(-(1i)*ph2))/(1-tc_low*a2*z_upp*exp(-(1i)*ph2));
    
    T_num = tw^2 + (a1^2)*(abs(z_low)^2)-2*a1*abs(z_low)*tw*cos(ph1-angle(z_low));
    T_denum = 1+(a1^2)*(abs(z_low)^2)*(tw^2)-2*a1*abs(z_low)*tw*cos(ph1-angle(z_low));
        
    T(counter1) = T_num/T_denum;
    counter1=counter1+1;
end

T_m = T_monos(tw, a, phi(1), phi(end), length(phi));
T_d = T_diplos(0.85, 0.8, a, a, phi(1), phi(end), length(phi), 0);

figure;
hold on;
% plot(phi, T, 'LineWidth', 1, 'Color', [0,1,0])
% plot(phi, T_m, 'LineWidth', 1, 'Color', [1,0,0])
% plot(phi, T_d, 'LineWidth', 1, 'Color', [0,0,1])

plot(phi, 10*log10(T), 'LineWidth', 1, 'Color', [0,1,0])
plot(phi, 10*log10(T_m), 'LineWidth', 1, 'Color', [1,0,0])
plot(phi, 10*log10(T_d), 'LineWidth', 1, 'Color', [0,0,1])

title('Συντελεστης Μεταδοσης vs φ')
xlabel('φ')
ylabel('T (dB)')
legend('T', 'T_m', 'T_d', 'Location', 'best');

hold off;


