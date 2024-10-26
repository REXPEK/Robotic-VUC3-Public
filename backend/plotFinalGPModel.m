function plotFinalGPModel(bo_results, var1, var2)    
    % Prepare data
    var1Data = bo_results.XTrace.(var1);
    var2Data = bo_results.XTrace.(var2);
    objData = bo_results.ObjectiveTrace;

    % Generate a scatter plot of the objective values
    figure;
    scatter3(var1Data, var2Data, objData, 36, objData, 'filled');
    colormap(jet);
    colorbar;
    xlabel(var1);
    ylabel(var2);
    zlabel('Objective Value');
    title(['Scatter Plot of Objective Values for ', var1, ' and ', var2]);
    grid on;
end
