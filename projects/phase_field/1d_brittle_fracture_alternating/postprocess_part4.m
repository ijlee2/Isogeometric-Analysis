%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine plots the errors in the strain, surface, and potential
%  energies, and finds the rate at which the errors decay to 0.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_postprocess_part4.m output
%  
%  Note that,
%  
%      order is the order of the phase field theory (2, 4, 6, 8)
%      path_to_assembly_directory is the path to the assembly directory
%      path_to_results_directory is the path to the results directory
%      path_to_outputs_directory is the path to the plots directory
%  
%  
%  Output:
%  
%  1. Convergence plots (.png files)
%--------------------------------------------------------------------------
function postprocess_part4()
    %----------------------------------------------------------------------
    %  Initialize the arrays
    %----------------------------------------------------------------------
    % Set the phase field theories
    orders = [2; 4; 6; 8];
    
    % Set the refinement levels
    numRefinements = [0; 1; 2; 3; 4];
    
    temp = zeros(size(numRefinements, 1), size(orders, 1));
    
    strain_energy_error    = temp;
    surface_energy_error   = temp;
    potential_energy_error = temp;
    
    clear temp;
    
    
    for j = 1 : size(orders, 1)
        % Load the global assembly file
        load(sprintf('../assembly_files/order%d/numRefinements%d/file_assembly_global', orders(j), numRefinements(1)), ...
             'material_L', ...
             'material_A', ...
             'material_G_c', ...
             'material_ell_0');
        
        % Load the BCs
        load(sprintf('../assembly_files/order%d/numRefinements%d/file_bc', orders(j), numRefinements(1)), 'u_L');
        
        
        %------------------------------------------------------------------
        %  Calculate the smallest and largest element sizes
        %------------------------------------------------------------------
        h_min = inf * ones(size(numRefinements, 1), 1);
        h_max = zeros(size(numRefinements, 1), 1);
        
        for i = 1 : size(numRefinements, 1)
            % Load the patch assembly file
            load(sprintf('../assembly_files/order%d/numRefinements%d/file_assembly_patch%d', orders(j), numRefinements(i), 1), 'elementSizes1');
            
            h_min(i) = min(h_min(i), min(elementSizes1));
            h_max(i) = max(h_max(i), max(elementSizes1));
        end
        
        
        %------------------------------------------------------------------
        %  Calculate the errors in the energies
        %------------------------------------------------------------------
        % Load the error file
        load(sprintf('../results/order%d/outputs/errors_and_energies', orders(j)), ...
             'strain_energy' , ...
             'surface_energy', ...
             'potential_energy');
        
        strain_energy_error(:, j)    = abs(strain_energy - 0);
        surface_energy_error(:, j)   = abs(surface_energy - material_G_c * material_A);
        potential_energy_error(:, j) = abs(potential_energy - (0 + material_G_c * material_A));
        
        
        % Plot the strain energy error
        figure(1);
        h_str(j) = loglog(h_min, strain_energy_error(:, j)   , '-s', 'LineWidth', 3, 'MarkerSize', 12); hold on;
        
        % Plot the surface energy error
        figure(2);
        h_sur(j) = loglog(h_min, surface_energy_error(:, j)  , '-s', 'LineWidth', 3, 'MarkerSize', 12); hold on;
        
        % Plot the potential energy error
        figure(3);
        h_pot(j) = loglog(h_min, potential_energy_error(:, j), '-s', 'LineWidth', 3, 'MarkerSize', 12); hold on;
    end
    
    
    path_to_plots_directory = '../results/plots/';
    if ~exist(path_to_plots_directory, 'dir')
        mkdir(path_to_plots_directory);
    end
    
    % Plot the strain energy error
    figure(1);
    set(h_str(1), 'Color', [0, 0, 0]);
    set(h_str(2), 'Color', [.4, .8, .4]);
    set(h_str(3), 'Color', 'b');
    set(h_str(4), 'Color', 'r');
    
    xlabel('h', 'FontSize', 90);
    ylabel('strain energy error', 'FontSize', 90);
    legend(h_str, {' 2nd order theory', ' 4th order theory', ' 6th order theory', ' 8th order theory'}, 'Location', 'Southeast');
    axis([1e-6, 1e-3, 1e-7, 1e0]);
    axis square;
    grid on;
    set(gca, 'FontSize', 54, 'YTick', [1e-7; 1e-6; 1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1e0]);
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf('%sstrain_energy_error.png', path_to_plots_directory));
    
    % Plot the surface energy error
    figure(2);
    set(h_sur(1), 'Color', [0, 0, 0]);
    set(h_sur(2), 'Color', [.4, .8, .4]);
    set(h_sur(3), 'Color', 'b');
    set(h_sur(4), 'Color', 'r');
    
    xlabel('h', 'FontSize', 90);
    ylabel('surface energy error', 'FontSize', 90);
    legend(h_sur, {' 2nd order theory', ' 4th order theory', ' 6th order theory', ' 8th order theory'}, 'Location', 'Southeast');
    axis([1e-6, 1e-3, 1e-7, 1e0]);
    axis square;
    grid on;
    set(gca, 'FontSize', 54, 'YTick', [1e-7; 1e-6; 1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1e0]);
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf('%ssurface_energy_error.png', path_to_plots_directory));
    
    % Plot the potential energy error
    figure(3);
    set(h_pot(1), 'Color', [0, 0, 0]);
    set(h_pot(2), 'Color', [.4, .8, .4]);
    set(h_pot(3), 'Color', 'b');
    set(h_pot(4), 'Color', 'r');
    
    xlabel('h', 'FontSize', 90);
    ylabel('potential energy error', 'FontSize', 90);
    legend(h_pot, {' 2nd order theory', ' 4th order theory', ' 6th order theory', ' 8th order theory'}, 'Location', 'Southeast');
    axis([1e-6, 1e-3, 1e-7, 1e0]);
    axis square;
    grid on;
    set(gca, 'FontSize', 54, 'YTick', [1e-7; 1e-6; 1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1e0]);
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf('%spotential_energy_error.png', path_to_plots_directory));
end
