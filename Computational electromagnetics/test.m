%contourf(real(E_gamma))
figure;
disp(Gamma)
El = -imag(E(:, end-1));
El = reshape(El, [], N);
p = max(max(abs(El)));
contourf(El);
hold on;
x_min = 37;
x_max = 73;
y_min = 37;
y_max = 73;

% Draw vertical lines (yline) from y_min to y_max at x = 37 and x = 73
line([37 37], [y_min y_max], 'color', 'r', 'LineWidth', 2)
line([73 73], [y_min y_max], 'color', 'r', 'LineWidth', 2)

% Draw horizontal lines (xline) from x_min to x_max at y = 37 and y = 73
yline(37, 'color', 'g', 'LineWidth', 2)
line([x_min x_max], [73 73], 'color', 'r', 'LineWidth', 2)
% line([0, grid_length], [y_pos1, y_pos1], 'Color', 'k', 'LineWidth', 2);
% line([0, grid_length], [y_pos2, y_pos2], 'Color', 'k', 'LineWidth', 2);
% 
% % Draw vertical lines
% line([x_pos1, x_pos1], [0, grid_length], 'Color', 'k', 'LineWidth', 2);
% line([x_pos2, x_pos2], [0, grid_length], 'Color', 'k', 'LineWidth', 2);

% 
% % Circle parameters
% radius = 30;
% 
% % Calculate the center of the grid
% center_x = N / 2;
% center_y = N / 2;
% 
% % Create the circle
% theta = linspace(0, 2*pi, 100);
% circle_x = center_x + radius * cos(theta);
% circle_y = center_y + radius * sin(theta);
% 
% % Plot the red circle
% plot(circle_x, circle_y, 'r', 'LineWidth', 2);
% 
% % Set axis labels and title
% xlabel('X axis');
% ylabel('Y axis');
% 
% % Ensure the aspect ratio is equal
% axis equal;
% 
% % Hold off to stop adding to the plot
% hold off;