%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine plots (1) the displacement field, (2) the derivatives of
%  the phase field, and (3) the stress field at the last Newton's iteration
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
function postprocess_part2(order, path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
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
                  0.40, 0.70, 0.90; ...
                  0.00, 0.00, 0.00];
    
    
    % Number of points that we consider on each element
    numPointsPerElement = 11;
    
    % Points on the Bernstein domain
    t = linspace(0, 1, numPointsPerElement)';
    
    % Some useful constants for evaluations
    p1 = order / 2;
    constant_p1p1 = p1 + 1;
    
    % Evaluate the Bernstein polynomials
    Bernstein1_der0 = zeros(constant_p1p1, numPointsPerElement);
    Bernstein1_der1 = zeros(constant_p1p1, numPointsPerElement);
    Bernstein1_der2 = zeros(constant_p1p1, numPointsPerElement);
    Bernstein1_der3 = zeros(constant_p1p1, numPointsPerElement);
    Bernstein1_der4 = zeros(constant_p1p1, numPointsPerElement);

    for q_e = 1 : numPointsPerElement
        temp = eval_1d_bernstein_der(t(q_e), p1);

        Bernstein1_der0(:, q_e) = temp(:, 1);
        Bernstein1_der1(:, q_e) = temp(:, 2);
        if (order >= 4)
            Bernstein1_der2(:, q_e) = temp(:, 3);
        end
        if (order >= 6)
            Bernstein1_der3(:, q_e) = temp(:, 4);
        end
        if (order >= 8)
            Bernstein1_der4(:, q_e) = temp(:, 5);
        end
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
             'maxNewtonsMethod');
        
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
        load(sprintf('%snumRefinements%d/file_results_iteration%d', path_to_results_directory, numRefinements(i), maxNewtonsMethod), 'u1', 'u2');
        
        
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
        u2_der0 = temp;
        u2_der1 = temp;
        if (order >= 4)
            u2_der2 = temp;
        end
        if (order >= 6)
            u2_der3 = temp;
        end
        if (order >= 8)
            u2_der4 = temp;
        end
        
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
            basis_der0     =             bezierExtractions1_e  * Bernstein1_der0;
            basis_der1     = (dt_dxi   * bezierExtractions1_e) * Bernstein1_der1;
            if (order >= 4)
                basis_der2 = (dt_dxi^2 * bezierExtractions1_e) * Bernstein1_der2;
            end
            if (order >= 6)
                basis_der3 = (dt_dxi^3 * bezierExtractions1_e) * Bernstein1_der3;
            end
            if (order >= 8)
                basis_der4 = (dt_dxi^4 * bezierExtractions1_e) * Bernstein1_der4;
            end
            
            
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
            u2_e = u2(LM_array2(:, e))';
            
            
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
                
                basis_physical_der0     =                        basis_der0(:, q_e);
                basis_physical_der1     = JacobianMatrix_inv   * basis_der1(:, q_e);
                if (order >= 4)
                    basis_physical_der2 = JacobianMatrix_inv^2 * basis_der2(:, q_e);
                end
                if (order >= 6)
                    basis_physical_der3 = JacobianMatrix_inv^3 * basis_der3(:, q_e);
                end
                if (order >= 8)
                    basis_physical_der4 = JacobianMatrix_inv^4 * basis_der4(:, q_e);
                end
                
                
                %----------------------------------------------------------
                %  Find the point on the physical domain
                %----------------------------------------------------------
                x(q) = nodes_e * basis_physical_der0;
                
                
                %----------------------------------------------------------
                %  Evaluate the displacement field
                %----------------------------------------------------------
                u1_der0(q) = u1_e * basis_physical_der0;
                u1_der1(q) = u1_e * basis_physical_der1;
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                u2_der0(q)     = u2_e * basis_physical_der0;
                u2_der1(q)     = u2_e * basis_physical_der1;
                if (order >= 4)
                    u2_der2(q) = u2_e * basis_physical_der2;
                end
                if (order >= 6)
                    u2_der3(q) = u2_e * basis_physical_der3;
                end
                if (order >= 8)
                    u2_der4(q) = u2_e * basis_physical_der4;
                end
                
                
                % Increment the counter
                q = q + 1;
            end
        end
        
        
        %------------------------------------------------------------------
        %  Evaluate the stress field
        %------------------------------------------------------------------
        degradation_function = u2_der0.^2;
        
        stress_xx = degradation_function .* (material_E * u1_der1);
        
        
        %------------------------------------------------------------------
        %  Plot the exact solutions
        %------------------------------------------------------------------
        if (i == size(numRefinements, 1))
            % Plot the exact displacement field
            figure(1);
            u1_der0_exact = function_u_exact(x, u_L);
            h(6) = plot(x, u1_der0_exact, '--', 'Color', line_color(6, :), 'LineWidth', 3); hold on;
            
            % Plot the exact phase field
            figure(2);
            u2_der0_exact = function_c_exact(x, material_ell_0, order);
            h(12) = plot(x, u2_der0_exact, '--', 'Color', line_color(6, :), 'LineWidth', 3); hold on;
        end
        
        
        %------------------------------------------------------------------
        %  Plot the displacement field
        %------------------------------------------------------------------
        figure(1);
        h(i) = plot(x, u1_der0, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        
        
        %------------------------------------------------------------------
        %  Plot the phase field
        %------------------------------------------------------------------
        figure(2);
        h(i + 6) = plot(x, u2_der0, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        
        figure(3);
        h(i + 12) = plot(x, u2_der1, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        
        if (order >= 4)
            figure(4);
            h(i + 17) = plot(x, u2_der2, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        end
        
        if (order >= 6)
            figure(5);
            h(i + 22) = plot(x, u2_der3, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        end
        
        if (order >= 8)
            figure(6);
            h(i + 27) = plot(x, u2_der4, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        end
        
        
        %------------------------------------------------------------------
        %  Plot the stress field
        %------------------------------------------------------------------
        figure(7);
        h(i + 32) = semilogy(x, stress_xx, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
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
    
    legend(h(1 : 6), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements', ' analytical solution'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, u_L, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_u_der%d_order%d.png', path_to_outputs_directory, 0, order)));
    
    
    % Phase field, 0th derivative
    figure(2);
    
    xlabel('x', 'FontSize', 90);
    ylabel('c(x)', 'FontSize', 90);
    axis([xmin xmax -0.05 1.05]);
    grid on;
    
    legend(h(7 : 12), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements', ' analytical solution'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, 1, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_c_der%d_order%d.png', path_to_outputs_directory, 0, order)));
    
    
    % Phase field, 1st derivative
    figure(3);
    
    xlabel('x', 'FontSize', 90);
    ylabel('c''(x)', 'FontSize', 90);
    axis([xmin xmax -inf inf]);
    grid on;
    
    legend(h(13 : 17), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_c_der%d_order%d.png', path_to_outputs_directory, 1, order)));
    
    
    % Phase field, 2nd derivative
    if (order >= 4)
        figure(4);
        
        xlabel('x', 'FontSize', 90);
        ylabel('c''''(x)', 'FontSize', 90);
        axis([xmin xmax -inf inf]);
        grid on;
        
        legend(h(18 : 22), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
        set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
        print('-dpng', '-r300', sprintf(sprintf('%splot_c_der%d_order%d.png', path_to_outputs_directory, 2, order)));
    end
    
    
    % Phase field, 3rd derivative
    if (order >= 6)
        figure(5);
        
        xlabel('x', 'FontSize', 90);
        ylabel('c''''''(x)', 'FontSize', 90);
        axis([xmin xmax -inf inf]);
        grid on;
        
        legend(h(23 : 27), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
        set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
        print('-dpng', '-r300', sprintf(sprintf('%splot_c_der%d_order%d.png', path_to_outputs_directory, 3, order)));
    end
    
    
    % Phase field, 4th derivative
    if (order >= 8)
        figure(6);
        
        xlabel('x', 'FontSize', 90);
        ylabel('c''''''''(x)', 'FontSize', 90);
        axis([xmin xmax -inf inf]);
        grid on;
        
        legend(h(28 : 32), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
        set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
        print('-dpng', '-r300', sprintf(sprintf('%splot_c_der%d_order%d.png', path_to_outputs_directory, 4, order)));
    end
    
    
    % Stress field
    figure(7);
    
    xlabel('x', 'FontSize', 90);
    ylabel('sigmaxx(x)', 'FontSize', 90);
    axis([xmin xmax 1e-10 1e6]);
    grid on;
    
    legend(h(33 : 37), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_sigmaxx_order%d.png', path_to_outputs_directory, order)));
    
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Plot the fields over the entire domain
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
             'maxNewtonsMethod');
        
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
        load(sprintf('%snumRefinements%d/file_results_iteration%d', path_to_results_directory, numRefinements(i), maxNewtonsMethod), 'u1', 'u2');
        
        
        %------------------------------------------------------------------
        %  Set parameters for plotting
        %------------------------------------------------------------------
        % Initialize the points on the physical domain
        temp = zeros(numPointsPerElement * numElements1, 1);
        
        x       = temp;
        u1_der0 = temp;
        u1_der1 = temp;
        u2_der0 = temp;
        
        clear temp;
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e = e1)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Counter for the point
        q = 1;
        
        
        for e = 1 : numElements1
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
            u2_e = u2(LM_array2(:, e))';
            
            
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
                
                
                %----------------------------------------------------------
                %  Evaluate the phase field
                %----------------------------------------------------------
                u2_der0(q) = u2_e * basis_physical_der0;
                
                
                % Increment the counter
                q = q + 1;
            end
        end
        
        
        %------------------------------------------------------------------
        %  Evaluate the stress field
        %------------------------------------------------------------------
        degradation_function = u2_der0.^2;
        
        stress_xx = degradation_function .* (material_E * u1_der1);
        
        
        %------------------------------------------------------------------
        %  Plot the displacement field
        %------------------------------------------------------------------
        figure(8);
        h(i + 37) = plot(x, u1_der0, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        
        
        %------------------------------------------------------------------
        %  Plot the phase field
        %------------------------------------------------------------------
        figure(9);
        h(i + 42) = plot(x, u2_der0, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
        
        
        %------------------------------------------------------------------
        %  Plot the stress field
        %------------------------------------------------------------------
        figure(10);
        h(i + 47) = semilogy(x, stress_xx, '-', 'Color', line_color(i, :), 'LineWidth', 3); hold on;
    end
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    xmin = -material_L;
    xmax = material_L;
    
    
    % Displacement field
    figure(8);
    
    xlabel('x', 'FontSize', 90);
    ylabel('u(x)', 'FontSize', 90);
    axis([xmin xmax -.05*u_L 1.05*u_L]);
    grid on;
    
    legend(h(38 : 42), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, u_L, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_u_order%d (entire domain).png', path_to_outputs_directory, order)));
    
    
    % Phase field
    figure(9);
    
    xlabel('x', 'FontSize', 90);
    ylabel('c(x)', 'FontSize', 90);
    axis([xmin xmax -.05 1.05]);
    grid on;
    
    legend(h(43 : 47), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(0, 1, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_c_order%d (entire domain).png', path_to_outputs_directory, order)));
    
    
    % Stress field
    figure(10);
    
    xlabel('x', 'FontSize', 90);
    ylabel('sigmaxx(x)', 'FontSize', 90);
    axis([xmin xmax 1e-10 1e6]);
    grid on;
    
    legend(h(48 : 52), {' 2000 elements', ' 4000 elements', ' 8000 elements', ' 16000 elements', ' 32000 elements'}, 'Location', 'Southeast');
    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf(sprintf('%splot_sigmaxx_order%d (entire domain).png', path_to_outputs_directory, order)));
end


function u = function_u_exact(x, u_L)
    u = u_L * (sign(x)/2 + 1/2);
end


function c = function_c_exact(x, ell_0, order)
    % Normalize the physical domain
    y = abs(x)/ell_0;
    
    switch order
        case 2
            c = 1 - exp(-1/2*y);
            
        case 4
            c = 1 - exp(-y) .* (1 + y);
            
        case 6
            c = 1 - exp(-3/2*y) .* (1 + 3/2*y + 9/8*y.^2);
            
        case 8
            c = 1 - exp(-2*y) .* (1 + 2*y + 2*y.^2 + 4/3*y.^3);
            
    end
end
