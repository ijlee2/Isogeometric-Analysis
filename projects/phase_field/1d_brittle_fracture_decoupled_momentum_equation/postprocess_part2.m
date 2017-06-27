%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine plots (1) the displacement field, (2) the derivatives of
%  the phase field, and (3) the stress field at the last alternation step
%  for each refinement level. We zoom in to see what happens near the crack.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_postprocess_part2.m output
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
%  1. Plots of the derivatives of the displacement and phase fields (.png files)
%--------------------------------------------------------------------------
function postprocess_part2(order, p, path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
    close all;
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    % Set the refinement levels
    numRefinements = [0; 1; 2; 3; 4];
    
    % Set the line colors for each refinement level
    line_color = [0.70, 0.70, 0.70; ...
                  0.70, 0.35, 0.45; ...
                  0.60, 0.55, 0.15; ...
                  0.30, 0.70, 0.50; ...
                  0.40, 0.70, 0.90];
    
    
    % Number of points that we consider on each element
    numPointsPerElement = 11;
    
    % Points on the Bernstein domain
    t = linspace(0, 1, numPointsPerElement)';
    
    % Some useful constants for evaluations
    p1 = p;
    constant_p1p1 = p1 + 1;
    
    % Evaluate the Bernstein polynomials
    Bernstein1_der0 = zeros(constant_p1p1, numPointsPerElement);
    Bernstein1_der1 = zeros(constant_p1p1, numPointsPerElement);

    for q_e = 1 : numPointsPerElement
        temp = eval_1d_bernstein_der(t(q_e), p1);

        Bernstein1_der0(:, q_e) = temp(:, 1);
        Bernstein1_der1(:, q_e) = temp(:, 2);
    end

    clear t temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Plot the fields
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for i = size(numRefinements, 1) : -1 : 1
        % Load the global assembly file
        load(sprintf('%snumRefinements%d/file_assembly_global', path_to_assembly_directory, numRefinements(i)), ...
             'material_L', ...
             'material_A', ...
             'material_G_c', ...
             'material_ell_0', ...
             ...
             'maxAlternations');
        
        % Load the patch assembly file
        load(sprintf('%snumRefinements%d/file_assembly_patch%d', path_to_assembly_directory, numRefinements(i), 1), ...
             'material_E'        , ...
             'nodes'             , ...
             'bezierExtractions1', ...
             'numElements1'      , ...
             'elementSizes1'     , ...
             'IEN_array');
        
        % Load the BCs
        load(sprintf('%snumRefinements%d/file_bc', path_to_assembly_directory, numRefinements(i)), 'LM_array1', 'LM_array2', 'u_L');
        
        % Load the results file
        load(sprintf('%snumRefinements%d/file_results_alternation%d', path_to_results_directory, numRefinements(i), 0));
        
        
        %------------------------------------------------------------------
        %  Set parameters for plotting
        %------------------------------------------------------------------
        % Specify the elements to zoom in on
        elementIndex_begin = numElements1 / 2 - (200 * 2^numRefinements(i));
        elementIndex_end   = numElements1 / 2 + (200 * 2^numRefinements(i));
        numElementsToConsider = elementIndex_end - elementIndex_begin + 1;
        
        % Initialize the points on the physical domain
        temp = zeros(numPointsPerElement * numElementsToConsider, 1);
        
        x       = temp;
        u1_der0 = temp;
        u1_der1 = temp;
        
        clear temp;
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Counter for the point
        q = 1;
        
        
        for e = elementIndex_begin : elementIndex_end
            % Find the Bezier extraction matrix
            bezierExtractions1_e = bezierExtractions1(:, :, e);
            
            % Evaluate the map derivative dt/dxi (constant)
            dt_dxi = 1 / elementSizes1(e);
            
            
            %--------------------------------------------------------------
            %  Evaluate the basis functions in the parametric domain
            %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            % Find the positions of the nodes
            nodes_e = nodes(IEN_array(:, e), :)';
            
            % Evaluate the B-splines
            basis_der0 =           bezierExtractions1_e  * Bernstein1_der0;
            basis_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
            
            
            %--------------------------------------------------------------
            %  Evaluate the map derivatives (x = x1)
            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
            %--------------------------------------------------------------
            dx_dxi = nodes_e * basis_der1;
            
            
            %--------------------------------------------------------------
            %  Loop over points
            %--------------------------------------------------------------
            % Find the coefficients of the solution vectors
            u1_e = u1(LM_array1(:, e))';
            
            
            for q_e = 1 : numPointsPerElement
                %----------------------------------------------------------
                %  Form the Jacobian matrix
                %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                %----------------------------------------------------------
                JacobianMatrix = dx_dxi(:, q_e);
                
                % Evaluate the Jacobian
                Jacobian = JacobianMatrix / dt_dxi;
                
                if (Jacobian <= 0)
                    fprintf('\n');
                    fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                    
                    quit;
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the basis functions in the physical domain
                %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                %----------------------------------------------------------
                JacobianMatrix_inv = 1 / JacobianMatrix;
                
                basis_physical_der0 =                      basis_der0(:, q_e);
                basis_physical_der1 = JacobianMatrix_inv * basis_der1(:, q_e);
                
                
                %----------------------------------------------------------
                %  Find the point on the physical domain
                %----------------------------------------------------------
                x(q) = nodes_e * basis_physical_der0;
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u1_der0(q) = u1_e * basis_physical_der0;
                u1_der1(q) = u1_e * basis_physical_der1;
                
                
                % Increment the counter
                q = q + 1;
            end
        end
        
        
        %------------------------------------------------------------------
        %  Evaluate the stress field
        %------------------------------------------------------------------
        y = abs(x) / material_ell_0;
        
        switch order
            case 2
                c_der0 = 1 - exp(-1/2*y);
                
            case 4
                c_der0 = 1 - exp(-y) .* (1 + y);
                
            case 6
                c_der0 = 1 - exp(-3/2*y) .* (1 + 3/2*y + 9/8*y.^2);
                
            case 8
                c_der0 = 1 - exp(-2*y) .* (1 + 2*y + 2*y.^2 + 4/3*y.^3);
                
        end
                
        
        stress_xx = c_der0.^2 .* (material_E * u1_der1);
        
        
        %------------------------------------------------------------------
        %  Plot the displacement field
        %------------------------------------------------------------------
        figure(1);
        h(i) = plot(x, u1_der0, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        
        
        %------------------------------------------------------------------
        %  Plot the stress field
        %------------------------------------------------------------------
        figure(2);
        h(i + 5) = semilogy(x, stress_xx, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Plot the fields
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    xmin = -1;
    xmax = 1;
    
    
    % Displacement field
    figure(1);
    
    xlabel('x', 'FontSize', 90);
    ylabel('u(x)', 'FontSize', 90);
    axis([xmin xmax -0.05*u_L 1.05*u_L]);
    grid on;
    
    legend(h(1 : 5), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, u_L, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_u_der0.png', path_to_outputs_directory)));
    
    
    % Stress field
    figure(2);
    
    xlabel('x', 'FontSize', 90);
    ylabel('sigmaxx(x)', 'FontSize', 90);
    axis([xmin xmax 1e-15 1e6]);
    grid on;
    
    legend(h(6 : 10), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_sigmaxx.png', path_to_outputs_directory)));
end
