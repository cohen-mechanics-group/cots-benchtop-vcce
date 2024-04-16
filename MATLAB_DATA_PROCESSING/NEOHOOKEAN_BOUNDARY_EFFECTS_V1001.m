% Brendan M Unikewicz, PhD Student
% Andre M Pincot, MSc Student
% Tal Cohen, Asc. Professor
% MIT, Dept. Mechanical Engineering
% MIT, Dept. Civil & Environmental Engineering
% Date of Creation: 03/26/2024
% Code Purpose: checking boundary effects for neo-Hookean material models

%% Code Start
clear; close all; clc; 

%% Define Colors
fontSize = 14; 
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], ...
    [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], ...
    [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330], [0.6350 0.0780 0.1840]};

%% Start Boundary Effects
% Define stretch range
lambda = linspace(0.05, 50, 10000);

% Define B/A ratios
BA_ratios = [5.0, 10.0, 15.0, 20.0, 1e6]; % Including infinity
BA_ratios_legend = {'5'; '10'; '15'; '20'; '$\infty$'};

% Calculate P/E for each B/A ratio
P_E = zeros(length(BA_ratios), length(lambda));
for i = 1:length(BA_ratios)
    if isinf(BA_ratios(i))
        % For B/A -> infinity, P/E is always 0
        P_E(i, :) = zeros(1, length(lambda));
    else
        % Calculate lambda_b using the provided equation
        A_over_B = 1/BA_ratios(i); % Since B/A is provided, we invert it
        lambda_b = (1 + (lambda.^3 - 1).*(A_over_B^3)).^(1/3);
        % Calculate P/E using the values of lambda_b and lambda
        P_E(i, :) = (1/6) .* (lambda_b.^-4 + 4.*lambda_b.^-1 - lambda.^-4 - 4.*lambda.^-1);
    end
    plot(lambda, P_E(i, :),'Color',colors{i},'LineWidth',1.5);
    ax = gca; ax.TickLabelInterpreter = 'latex'; ax.XAxis.TickLabelInterpreter = 'latex'; ax.YAxis.TickLabelInterpreter = 'latex'; 
    hold on; % This will keep the current plot so the next one doesn't overwrite it
end
asymptote = repmat((5/6),1,length(lambda));
plot(lambda,asymptote,'--k','LineWidth',1.5);

txt_AL = 'Asymptotic Limit: 5/6';
text(4.75,0.916,txt_AL,'HorizontalAlignment','right','Interpreter','latex','FontSize',13);

xlabel('Stretch ($\lambda$)','Interpreter','latex','FontSize',14); xlim([1 5]);
ylabel('p/E','Interpreter','latex','FontSize',14); ylim([0 1]);
legend(['B/A = ' BA_ratios_legend{1}],['B/A = ' BA_ratios_legend{2}],['B/A = ' BA_ratios_legend{3}],['B/A = ' BA_ratios_legend{4}],['B/A $\rightarrow$ ' BA_ratios_legend{5}],'Interpreter','latex','Location','southeast','Fontsize',13);
hold off;
