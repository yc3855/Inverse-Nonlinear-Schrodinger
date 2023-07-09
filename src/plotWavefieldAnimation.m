function plotWavefieldAnimation(u_real, u_imag, step, option)
    % Get the size of the data
    [Nx, Ny, Nt] = size(u_real);
    
    % Create a figure
    figure;
    
    % Loop through the frames and plot
    for t = 1:step:Nt
        % Get the current frame of the wavefield
        u_t = u_real(:, :, t) + 1i * u_imag(:, :, t);
        
        % Determine the plot title based on the option
        if strcmp(option, 'amplitude')
            titleStr = 'Amplitude';
            plotData = abs(u_t);
        elseif strcmp(option, 'phase')
            titleStr = 'Phase';
            plotData = angle(u_t);
        else
            error('Invalid option. Please choose either "amplitude" or "phase"');
        end
        
        % Plot the contour
        contourf(plotData, 'LineWidth', 1.5);
        title(sprintf('%s (t = %d)', titleStr, t));
        colorbar;
        axis equal;
        pause(0.1);  % Pause between frames
        
        if t ~= Nt
            clf;  % Clear the figure for the next frame
        end
    end
end
