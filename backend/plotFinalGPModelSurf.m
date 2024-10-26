function plotFinalGPModelSurf(bo_results, var1, var2)    
    % Prepare data
    var1Data = bo_results.XTrace.(var1);
    var2Data = bo_results.XTrace.(var2);
    objData = bo_results.ObjectiveTrace;

    % Define the grid for var1 and var2
    [var1Grid, var2Grid] = meshgrid(linspace(min(var1Data), max(var1Data), 100), ...
                                    linspace(min(var2Data), max(var2Data), 100));
    
    % Interpolate objective data over the grid
    objGrid = griddata(var1Data, var2Data, objData, var1Grid, var2Grid, 'cubic');

    % Generate a surface plot of the objective values
    figure;
    surf(var1Grid, var2Grid, objGrid);
    colormap(jet);
    colorbar;
    xlabel(var1);
    ylabel(var2);
    zlabel('Objective Value');
    title(['Surface Plot of Objective Values for ', var1, ' and ', var2]);
    grid on;
end
