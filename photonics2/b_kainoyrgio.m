a = 0.85;
a1 = a;
a2 = a;
phi = pi:0.01*pi:5*pi;
tw = 0.85;
tc_values = [0.8, 0.9, 0.95];

counter2=1;

T = zeros(length(tc_values), length(phi));

for vv=tc_values    
    T(counter2, :) = T_diplos(tw, vv, a, a, phi(1), phi(end), length(phi), 0);
    counter2= counter2+1;
end

%gia to memonomeno sintonisti
T_a = zeros(1, length(phi));

for i=phi
    T_a(:) = T_monos(tw, a, phi(1), phi(end), length(phi));
end



figure;
hold on;

% plot(phi, T(1, :), 'LineWidth', 1, 'Color', [0,1,0])
% plot(phi, T(2, :), 'LineWidth', 1, 'Color',  [1,0,0])
% plot(phi, T(3, :), 'LineWidth', 1, 'Color', [0,0, 1])
% plot(phi, T_a, '--', 'LineWidth', 1, 'Color', [0,0, 0])

plot(phi, 10*log10(T(1, :)), 'LineWidth', 1, 'Color', [0,1,0])
plot(phi, 10*log10(T(2, :)), 'LineWidth', 1, 'Color',  [1,0,0])
plot(phi, 10*log10(T(3, :)), 'LineWidth', 1, 'Color', [0,0, 1])
plot(phi, 10*log10(T_a), '--', 'LineWidth', 1, 'Color', [0,0, 0])

title('Συντελεστης Μεταδοσης vs φ')
xlabel('φ')
ylabel('T(dB)')
legend('t_c = 0.8', 't_c = 0.9', 't_c = 0.95', 'T_m', 'Location', 'best');

hold off;
