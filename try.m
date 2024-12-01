% Data for the top plot
x1 = 0:7;
y1 = [0, 2, 4, 0, -2, 0, 4, -2];

% Data for the bottom plot
x2 = 0:7; 
y2 = [2, -4, 4, -2, 0, 0, 0, 0]; 

% Create a figure
figure;

% Top plot
subplot(2, 1, 1);
plot(x1, y1, 'r', 'LineWidth', 2); % Red line
grid on;
set(gca, 'Color', 'w', 'GridColor', 'g', 'GridLineStyle', '-', 'GridAlpha', 0.5);
set(gca, 'XColor', 'k', 'YColor', 'k'); % Black axes labels
xlim([0 7]);
ylim([-4 6]);

% Bottom plot
subplot(2, 1, 2);
plot(x2, y2, 'b', 'LineWidth', 2); % Blue line
grid on;
set(gca, 'Color', 'w', 'GridColor', 'g', 'GridLineStyle', '-', 'GridAlpha', 0.5);
set(gca, 'XColor', 'k', 'YColor', 'k'); % Black axes labels
xlim([0 7]);
ylim([-4 4]);

% Adjust figure background color to white
set(gcf, 'Color', 'w');
