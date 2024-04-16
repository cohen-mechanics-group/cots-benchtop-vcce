function [fitParams, fullFit_region, fullFit_dataY, fitParams_optimized] = neoHookean_fitting(percentile, pointPairs, eGuess, dataX, dataY, dataX_fullFit,total_aBar)

    % Define your model function
    model = @(C, x) C(2) .* ((5/6) - (2./(3.*(x./C(1)))) - (1./(6.*(x./C(1)).^4)));

    numCells = size(pointPairs, 2);
    fitParams = cell(1, numCells);
    fullFit_region = cell(1, numCells);
    
    for f = 1:numCells
        % Store results for pointPairData
        fitParams_pointPair = struct();
        fullFit_pointPair = struct();
        
        fieldNames = fieldnames(pointPairs{f});

        % Get the 90th percentile slope for current cell
        allSlopes = arrayfun(@(x) pointPairs{f}.(x{1}).slope, fieldNames);
        thresholdSlope = prctile(allSlopes, percentile);
        
        for idx = 1:length(fieldNames)
            fieldName = fieldNames{idx};
            currentData = pointPairs{f}.(fieldName);
            
            % Continue to next iteration if slope is below threshold
            if currentData.slope < thresholdSlope
                continue;
            end
            
            disp(['Type of currentData for ', fieldName, ': ', class(currentData)]);

            initial_guess = [currentData.x_intercept, eGuess];
            
            bound_defect = 10; 
            lb = [currentData.x_intercept*(-bound_defect), 1];
            ub = [currentData.x_intercept*(bound_defect), eGuess*10];
    
            fitParams_optimized = lsqcurvefit(model, initial_guess, dataX{f}/currentData.x_intercept, dataY{f}, lb, ub);
            disp('Running!')
            fitParams_pointPair.(fieldName) = fitParams_optimized;
            
            % Store the y_theory for the whole dataset
            fullFit_pointPair.(fieldName) = model(fitParams_optimized, dataX_fullFit{f});

            % Define your model function
            fullFit_dataY{f} = fitParams_optimized(2) .* ((5/6) - (2./(3.*((total_aBar{f}/currentData.x_intercept)./fitParams_optimized(1)))) - (1./(6.*((total_aBar{f}/currentData.x_intercept)./fitParams_optimized(1)).^4)));
        end
    
        % Store results in cell array
        fitParams{f} = struct('pointPairData', fitParams_pointPair);
        fullFit_region{f} = struct('pointPairData', fullFit_pointPair);

        % Define your model function
        % fullFit_dataY{f} = fitParams_optimized(2) .* ((5/6) - (2./(3.*(total_aBar{f}./fitParams_optimized(1)))) - (1./(6.*(total_aBar{f}./fitParams_optimized(1)).^4)));
    end
end



