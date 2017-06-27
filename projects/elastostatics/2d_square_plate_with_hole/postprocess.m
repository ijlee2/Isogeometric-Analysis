%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine plots the 2D specimen after deformation.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_postprocess.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%      path_to_outputs_directory is the path to the plots directory
%  
%  
%  Output:
%  
%  1. Plot of the displacement fields (.png file)
%--------------------------------------------------------------------------
function postprocess(path_to_assembly_directory, path_to_results_directory, path_to_outputs_directory)
    close all;
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches', ...
         'numNodesBeforePatch', ...
         'numDOFsPerNode', ...
         'GN_array');
    
    % Load the BC file
    load(sprintf('%sfile_bc', path_to_assembly_directory), ...
         'ID_array');
    
    % Load the results file
    load(sprintf('%sfile_results', path_to_results_directory), ...
         'u');
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over patches
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        % Load the patch assembly file
        load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, p), ...
             'elementType', ...
             'p1', ...
             'p2', ...
             'nodes', ...
             'numNodesPerElement', ...
             'bezierExtractions1', ...
             'bezierExtractions2', ...
             'numElements1', ...
             'numElements2', ...
             'IEN_array');
        
        
        switch (lower(elementType))
            case {'bspline', 'b-spline'}
                elementType = 0;
                
            case 'nurbs'
                elementType = 1;
                
        end
        
        % Some useful constants
        constant_p1p1 = p1 + 1;
        constant_p2p1 = p2 + 1;
        
        
        %------------------------------------------------------------------
        %  Build the LM array
        %------------------------------------------------------------------
        LM_array = build_lm_array(IEN_array, ID_array, GN_array, numNodesBeforePatch(p));
        
        
        %------------------------------------------------------------------
        %  Set parameters for plotting
        %------------------------------------------------------------------
        % Number of points where the fields are evaluated
        numPointsPerElementCoordinate = 9;
        
        % Points on the Bernstein domain
        t = linspace(0, 1, numPointsPerElementCoordinate)';
        
        % Evaluate the Bernstein polynomials
        Bernstein1_der0 = zeros(constant_p1p1, numPointsPerElementCoordinate);
        Bernstein2_der0 = zeros(constant_p2p1, numPointsPerElementCoordinate);
        
        for j = 1 : numPointsPerElementCoordinate
            Bernstein1_der0(:, j) = eval_1d_bernstein(t(j), p1);
            Bernstein2_der0(:, j) = eval_1d_bernstein(t(j), p2);
        end
        
        % Initialize the points on the physical domain
        x = zeros(numPointsPerElementCoordinate, numPointsPerElementCoordinate);
        y = zeros(numPointsPerElementCoordinate, numPointsPerElementCoordinate);
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Counter for the element
        e = 1;
        
        for e2 = 1 : numElements2
            % Find the Bezier extraction matrix
            bezierExtractions2_e = bezierExtractions2(:, :, e2);
            
            % Evaluate the univariate B-splines
            Bspline2_der0 = bezierExtractions2_e * Bernstein2_der0;
            
            

            for e1 = 1 : numElements1
                % Find the Bezier extraction matrix
                bezierExtractions1_e = bezierExtractions1(:, :, e1);
                
                % Evaluate the univariate B-splines
                Bspline1_der0 = bezierExtractions1_e * Bernstein1_der0;
                
                
                %----------------------------------------------------------
                %  Evaluate the B-splines
                %----------------------------------------------------------
                Bspline_der00 = kron(Bspline2_der0, Bspline1_der0);
                
                
                %----------------------------------------------------------
                %  Find the positions of the nodes after deformation
                %----------------------------------------------------------
                % Find the positions of the nodes
                nodes_e = nodes(IEN_array(:, e), :)';
                
                % Find the displacements of the nodes
                u_e = u(LM_array(:, e));
                
                % Counter for the displacements of the nodes
                index = 1;
                
                for j = 1 : numNodesPerElement
                    for i = 1 : numDOFsPerNode
                        if (elementType == 0)
                            nodes_e(i, j) = nodes_e(i, j) + u_e(index);
                        
                        else
                            % Project the NURBS to the B-spline in the
                            % 3-dimensional space
                            nodes_e(i, j) = (nodes_e(i, j) + u_e(index)) * nodes_e(3, j);
                            
                        end
                        
                        index = index + 1;
                    end
                end
                
                
                %----------------------------------------------------------
                %  Loop over the points
                %----------------------------------------------------------
                % Counter for the point
                j = 1;
                
                for j2 = 1 : numPointsPerElementCoordinate
                    for j1 = 1 : numPointsPerElementCoordinate
                        % Evaluate the basis functions in the physical domain
                        basis_der00 = Bspline_der00(:, j);
                        
                        % Find the point on the physical domain
                        if (elementType == 0)
                            x(j1, j2) = nodes_e(1, :) * basis_der00;
                            y(j1, j2) = nodes_e(2, :) * basis_der00;
                            
                        elseif (elementType == 1)
                            denominator = nodes_e(3, :) * basis_der00;
                            
                            x(j1, j2) = (nodes_e(1, :) * basis_der00) / denominator;
                            y(j1, j2) = (nodes_e(2, :) * basis_der00) / denominator;
                            
                        end
                         
                        j = j + 1;
                    end
                end
                
                
                %----------------------------------------------------------
                %  Draw the B-spline surface
                %----------------------------------------------------------
                draw_2d_element(x, y, numPointsPerElementCoordinate);
                
                
                e = e + 1;
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over patches
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    xlabel('x', 'FontSize', 90);
    ylabel('y', 'FontSize', 90, 'Rotation', 0);
    axis([0 0.12 0 0.12]);
    axis square;
    grid on;
    set(gca, 'FontSize', 54);
%    set(gca, 'FontSize', 54, 'XTick', linspace(xmin, xmax, 11)', 'YTick', linspace(ymin, ymax, 11)');
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 40, 22.5]);
    print('-dpng', '-r300', sprintf('%splot_displacement.png', path_to_outputs_directory));
end


%--------------------------------------------------------------------------
%  This routine draws a contour of the 2D B-spline element. Points are
%  added throughout the element as a visual aid. Note that we draw these
%  points all at once for computational efficiency.
%--------------------------------------------------------------------------
function draw_2d_element(x, y, numPointsPerElementCoordinate)
    %----------------------------------------------------------------------
    %  Draw the points
    %----------------------------------------------------------------------
    % Number of points on an edge, excluding the endpoints
    numInteriorPointsPerEdge = (numPointsPerElementCoordinate - 1) / 4 - 1;
    
    % Initialize the points arrays
    points_on_edge = zeros(4 * numInteriorPointsPerEdge + 4, 2);
    points_in_element = zeros(numInteriorPointsPerEdge^2, 2);
    
    % Counters for the points arrays
    count_points_on_edge = 1;
    count_points_in_element = 1;
    
    % Select every other fourth point
    for j2 = 1 : 4 : numPointsPerElementCoordinate
        for j1 = 1 : 4 : numPointsPerElementCoordinate
            % Flag to check where the point exists in the volume
            flag = (j1 == 1 || j1 == numPointsPerElementCoordinate) + (j2 == 1 || j2 == numPointsPerElementCoordinate);
            
            % Mark the points on the edges of the outer surface
            if (flag >= 1)
                points_on_edge(count_points_on_edge, :) = [x(j1, j2), y(j1, j2)];
                count_points_on_edge = count_points_on_edge + 1;
                
            % Mark the points inside the area
            else
                points_in_element(count_points_in_element, :) = [x(j1, j2), y(j1, j2)];
                count_points_in_element = count_points_in_element + 1;
                
            end
        end
    end
    
    % Draw the points
    plot(points_on_edge(:, 1), points_on_edge(:, 2), ...
         's', 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.20; 0.50; 0.80], 'MarkerSize', 8); hold on;
    
    plot(points_in_element(:, 1), points_in_element(:, 2), ...
         'h', 'MarkerEdgeColor', [0.30; 0.90; 0.25], 'MarkerFaceColor', [0.30; 0.90; 0.25], 'MarkerSize', 6); hold on;
    
    clear points_on_edge points_in_element;
    
    
    %----------------------------------------------------------------------
    %  Draw the surfaces
    %----------------------------------------------------------------------
    % Draw the left surface
    j1 = 1;
    plot(x(j1, :), y(j1, :), ...
         '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
    
    % Draw the right surface
    j1 = numPointsPerElementCoordinate;
    plot(x(j1, :), y(j1, :), ...
         '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
    
    % Draw the bottom surface
    j2 = 1;
    plot(x(:, j2), y(:, j2), ...
         '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
    
    % Draw the top surface
    j2 = numPointsPerElementCoordinate;
    plot(x(:, j2), y(:, j2), ...
         '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
end