function [maxSlopeData, pointPairData, a_bar, kPa_interest, a_bar_interest, a_bar_split, indices_a_bar_interest, slopeData_experimental_a_bar,volData_experimental_totalInterest,slopeData_experimental_P,pData_experimental_totalInterest] = pointPairs_initialDefect(N_LOC_INPUT, V, kPa)

nL_to_mm3 = 0.001; 

% Needle Specifications
in_2_mm = 25.4; 
OD_25g_in = 0.020; 
ID_25g_in = 0.012; 

OD_25g_mm = OD_25g_in * in_2_mm;
OR_25g_mm = OD_25g_mm / 2; 
ID_25g_mm = ID_25g_in * in_2_mm;
IR_25g_mm = ID_25g_mm / 2; 

N_LOC = N_LOC_INPUT + 1; 

for f = 1:size(kPa, 1)
    %% Convert to 
    a_bar{f} = (((3/4) * (V{f} * nL_to_mm3)) / pi).^(1/3); % radius, mm

    %% FindPks & Region of Interest
    [pks{f},locs{f}] = findpeaks(kPa{f}); % all peaks
    [maxPeakValue{f}, index{f}] = max(pks{f});

    % equivalent system
    a_bar_split{f} = OR_25g_mm:((a_bar{f}(locs{f}(index{f})) - OR_25g_mm)/N_LOC):a_bar{f}(locs{f}(index{f})); 

    % The a_bar values you're interested in
    x_values{f} = [a_bar_split{f}]; 

    % Preallocate arrays for closest a_bar and corresponding kPa values
    indices_a_bar_interest{f} = zeros(size(x_values{f}));  %<--- New line here

    % Loop through the x_values
    for i = 1:length(x_values{f})
        % Compute absolute differences
        diff{f} = abs(a_bar{f} - x_values{f}(i));

        % Find index of minimum difference
        [~, idx1{f}] = min(diff{f});

        % Store closest a_bar and corresponding kPa
        a_bar_interest{f}(i) = a_bar{f}(idx1{f});
        kPa_interest{f}(i) = kPa{f}(idx1{f});
        indices_a_bar_interest{f}(i) = idx1{f}; %<--- New line here
    end 

    % Display results
    for i = 1:length(x_values{f})
        fprintf('Closest to %.2f: %.2f with corresponding kPa value: %.2f\n', ...
            x_values{f}(i), a_bar_interest{f}(i), kPa_interest{f}(i));
    end

    %% Lines of best fit over 2-pair combinations
    % Total number of points of interest
    N{f} = length(a_bar_interest{f});

% %     % Prepare a figure for the plots
% %     figure; hold on;

    % Initialize the struct to save pointPairData
    pointPairData{f} = struct();

    % Initialize the maximum slope and corresponding point pair
    maxSlope{f} = -inf;
    maxPointPair{f} = '';

    % Iterate over all combinations of two points
    for i = 1:N{f}
        for j = i+1:N{f}
            % Find indices where a_bar is between the two points of interest
            idx2 = (a_bar{f} >= a_bar_interest{f}(i) & a_bar{f} <= a_bar_interest{f}(j));

            % Extract the subset of the data between these points
            a_subset{f} = a_bar{f}(idx2);
            kpa_subset{f} = kPa{f}(idx2);

            % Fit a line to this subset of data
            p{f} = polyfit(a_subset{f}, kpa_subset{f}, 1);

            % Generate a line of best fit based on this fitted line
            yfit{f} = polyval(p{f}, a_subset{f});

% %             % Plot this line of best fit and label it
% %             subplot(2,1,1); 
% %             plot(a_bar{f},kPa{f},'k','LineWidth',1.3); hold on; 
% %             plot(a_subset{f}, yfit{f},'r','LineWidth',1.1); 
% %             xlabel('a'); ylabel('kPa'); grid on; 
% %             legend('P-a Curve','pointPairData', 'Location', 'best');
% %             % pointPair legendEntries are kind of wonky for legend()
            legendEntry = sprintf('pointPair%d%d', i, j);
% %             %legend(legendEntry);

            % Calculate the x-intercept for the line of best fit
            x_intercept{f} = -p{f}(2)/p{f}(1);

            % Save the slope, abscissa, ordinate, and x_intercept in a struct
            pointPairData{f}.(legendEntry).slope = p{f}(1);
            pointPairData{f}.(legendEntry).abscissa = a_subset{f};
            pointPairData{f}.(legendEntry).ordinate = yfit{f};
            pointPairData{f}.(legendEntry).x_intercept = x_intercept{f};
            
            % Update the maximum slope and corresponding point pair if the
            % current slope is greater than the maximum slope
            if p{f}(1) > maxSlope{f}
                maxSlope{f} = p{f}(1);
                maxPointPair{f} = legendEntry;
            end
        end
    end

    % 1. Extract the subset of data between indices_a_bar_interest(1) and indices_a_bar_interest(end)
    subset_P = kPa{f}(indices_a_bar_interest{f}(1):indices_a_bar_interest{f}(end));
    subset_a_bar = a_bar{f}(indices_a_bar_interest{f}(1):indices_a_bar_interest{f}(end));
    
    % 2. Calculate the middle index of this subset
    % middle_idx = round(length(subset_P) / 2);
    middle_idx = round(length(subset_P) * 0.55);
    
    % 3. Calculate 25% of the length of this subset
    quarter_length_low = round(0.2 * length(subset_P));
    quarter_length_high = round(0.3 * length(subset_P));
    
    % Calculate start and end indices for the 50% range in the subset
    start_index = max(1, middle_idx - quarter_length_low);
    end_index = min(length(subset_P), middle_idx + quarter_length_high);
    
    % 4. Extract the 50% data from the subset based on calculated indices
    slopeData_experimental_a_bar{f} = subset_a_bar(start_index:end_index);
    slopeData_experimental_P{f} = subset_P(start_index:end_index);

    volData_experimental_totalInterest{f} = a_bar{f}(indices_a_bar_interest{f}(1):indices_a_bar_interest{f}(end)); 
    pData_experimental_totalInterest{f} = kPa{f}(indices_a_bar_interest{f}(1):indices_a_bar_interest{f}(end)); 

    % Plot the line with the maximum slope
    maxSlopeData{f} = pointPairData{f}.(maxPointPair{f});

% %     subplot(2,1,2); 
% %     plot(maxSlopeData{f}.abscissa, maxSlopeData{f}.ordinate, 'LineWidth', 3, 'Color', 'm'); hold on; 
% % 
% %     % Display the point pair with the maximum slope
% %     fprintf('Point pair with the maximum slope: %s, slope = %.2f\n', maxPointPair{f}, maxSlope{f});
% % 
% %     % Plot the original data and points of interest
% %     plot(a_bar_interest{f}, kPa_interest{f}, 'or', 'LineWidth', 1.3);
% %     plot(a_bar{f}, kPa{f}, 'k', 'LineWidth', 1.3);
% %     xlabel('a'); ylabel('kPa');
% %     legend('maxSlope-pointPairSet','locations-pointsPairSets','P-a Curve','Location', 'best');
% %     grid on; 
% %     
% %     hold off; 
end

