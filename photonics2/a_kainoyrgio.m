tw = [0.6, 0.85, 0.95];
a = 0.85;
phi = pi:0.01*pi:5*pi;
counter1=1;
counter2=1;
T = zeros(length(tw), length(phi));

for tj=tw
    T(counter2, :) = T_monos(tj, a, phi(1), phi(end), length(phi));
    counter2=counter2+1;
end

figure;
hold on;
% plot(phi, T(1, :), 'LineWidth', 2, 'Color', [0,1,0])
% plot(phi, T(2, :), 'LineWidth', 2, 'Color',  [1,0,0])
% plot(phi, T(3, :), 'LineWidth', 2, 'Color', [0,0, 1])

plot(phi, 10*log10(T(1, :)), 'LineWidth', 1, 'Color', [0,1,0])
plot(phi, 10*log10(T(2, :)), 'LineWidth', 1, 'Color',  [1,0,0])
plot(phi, 10*log10(T(3, :)), 'LineWidth', 1, 'Color', [0,0, 1])

title('Συντελεστης Μεταδοσης vs φ')
xlabel('φ')
ylabel('T (dB)')
legend('t_w = 0.6', 't_w = 0.85', 't_w = 0.95', 'Location', 'best');

hold off;