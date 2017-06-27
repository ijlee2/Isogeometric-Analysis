%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine draws the NURBS surface that is specified by the knot
%  vectors of degrees p1, p2 and the nodes (control points) in 2D or 3D.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open. The nodes are assumed to be
%  ordered along the direction of 1, then along that of 2.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      draw_nurbs_surface(knots1, knots2, nodes, p1, p2);
%  
%  where,
%  
%      knots1, knots2 are the knot vectors (column vectors)
%      nodes is an array of nodes
%      p1, p2 are the degrees of the NURBS
%  
%  
%  Output:
%  
%  1. Plot of the NURBS surface
%--------------------------------------------------------------------------
function draw_nurbs_surface(knots1, knots2, nodes, p1, p2)
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    
    numKnots1 = size(knots1, 1);
    numKnots2 = size(knots2, 1);
    numNodes1 = numKnots1 - constant_p1p1;
    numNodes2 = numKnots2 - constant_p2p1;
    numNodes = numNodes1 * numNodes2;
    numDimensions = size(nodes, 2) - 1;
    
    % Check for errors
    if (numNodes ~= size(nodes, 1))
        fprintf('Error: For the given knot vectors, there should be %d nodes.\n\n', numNodes);
        
        return;
        
    elseif (numDimensions < 2 || numDimensions > 3)
        fprintf('Error: The nodes must be in the 2D or 3D space.\n\n');
        
        return;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Build the Bezier extraction matrices
    %----------------------------------------------------------------------
    [bezierExtractions1, nodeIndexShifts1, numElements1] = build_bezier_extraction(knots1, p1);
    [bezierExtractions2, nodeIndexShifts2, numElements2] = build_bezier_extraction(knots2, p2);
    
    
    %----------------------------------------------------------------------
    %  Build the IEN array
    %----------------------------------------------------------------------
    IEN_array = build_ien_array(knots1, knots2, [], p1, p2, [], nodeIndexShifts1, nodeIndexShifts2, []);
    
    
    %----------------------------------------------------------------------
    %  Project NURBS to the B-spline in the higher dimension
    %----------------------------------------------------------------------
    nodes_projectUp = project_up(nodes);
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    figure;
    
    % Number of points where the fields are evaluated
    if (min([numElements1; numElements2]) > 8)
        numPointsPerElementCoordinate = 9;
    else
        numPointsPerElementCoordinate = 17;
    end
    
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
    z = zeros(numPointsPerElementCoordinate, numPointsPerElementCoordinate);
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    % Counter for the element
    e = 1;
    
    for e2 = 1 : numElements2
        % Find the Bezier extraction matrix for direction 2
        bezierExtractions2_e = bezierExtractions2(:, :, e2);
        
        % Evaluate the B-splines for direction 2
        Bspline2_der0 = bezierExtractions2_e * Bernstein2_der0;
        
        
        for e1 = 1 : numElements1
            % Find the positions of the nodes
            nodes_e = nodes_projectUp(IEN_array(:, e), :)';
            
            % Find the Bezier extraction matrix for direction 1
            bezierExtractions1_e = bezierExtractions1(:, :, e1);
            
            % Evaluate the B-splines for direction 1
            Bspline1_der0 = bezierExtractions1_e * Bernstein1_der0;
            
            % Take the tensor product to get the B-splines in 2D
            Bspline_der0 = kron(Bspline2_der0, Bspline1_der0);
            
            
            %--------------------------------------------------------------
            %  Loop over the points
            %--------------------------------------------------------------
            % Counter for the point
            j = 1;
            
            for j2 = 1 : numPointsPerElementCoordinate
                for j1 = 1 : numPointsPerElementCoordinate
                    % Evaluate the basis functions for the physical domain
                    basis_der0 = Bspline_der0(:, j);
                    
                    % Find the weight
                    w = nodes_e(numDimensions + 1, :) * basis_der0;
                    
                    % Find the corresponding point on the physical domain
                    x(j1, j2) = (nodes_e(1, :) * basis_der0) / w;
                    
                    y(j1, j2) = (nodes_e(2, :) * basis_der0) / w;
                    
                    if (numDimensions == 3)
                        z(j1, j2) = (nodes_e(3, :) * basis_der0) / w;
                    end
                    
                    j = j + 1;
                end
            end
            
            
            %--------------------------------------------------------------
            %  Draw the NURBS surface
            %--------------------------------------------------------------
            %{
            % Uncomment this block to draw the control mesh
            nodes_e = nodes(IEN_array(:, e), :)';
            
            if (numDimensions == 2)
                mesh(reshape(nodes_e(1, :), [constant_p1p1, constant_p2p1]), ...
                     reshape(nodes_e(2, :), [constant_p1p1, constant_p2p1]), ...
                     zeros(constant_p1p1, constant_p2p1), ...
                     'EdgeAlpha', 0.5, 'EdgeColor', [0.75; 0.75; 0.75], 'FaceAlpha', 0.4, 'LineWidth', 2); hold on;
            else
                mesh(reshape(nodes_e(1, :), [constant_p1p1, constant_p2p1]), ...
                     reshape(nodes_e(2, :), [constant_p1p1, constant_p2p1]), ...
                     reshape(nodes_e(3, :), [constant_p1p1, constant_p2p1]), ...
                     'EdgeAlpha', 0.5, 'EdgeColor', [0.75; 0.75; 0.75], 'FaceAlpha', 0.4, 'LineWidth', 2); hold on;
            end
            %}
            
            draw_2d_element(x, y, z, numPointsPerElementCoordinate);
            
            
            e = e + 1;
        end
        
    end
    
    xlabel('x', 'FontSize', 90);
    ylabel('y', 'FontSize', 90);
    zlabel('z', 'FontSize', 90, 'Rotation', 0);
    axis square;
    grid on;
    set(gca, 'FontSize', 54);
    if (numDimensions == 2)
        ylabel('y', 'FontSize', 90, 'Rotation', 0);
        view(2);
    end
end


%--------------------------------------------------------------------------
%  This routine draws a contour of the 2D NURBS element. Points are
%  added throughout the element as a visual aid. Note that we draw these
%  points all at once for computational efficiency.
%--------------------------------------------------------------------------
function draw_2d_element(x, y, z, numPointsPerElementCoordinate)
    %----------------------------------------------------------------------
    %  Draw the points
    %----------------------------------------------------------------------
    % Number of points on an edge, excluding the endpoints
    numInteriorPointsPerEdge = (numPointsPerElementCoordinate - 1) / 4 - 1;
    
    % Initialize the points arrays
    points_on_edge = zeros(4 * numInteriorPointsPerEdge + 4, 3);
    points_in_element = zeros(numInteriorPointsPerEdge^2, 3);
    
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
                points_on_edge(count_points_on_edge, :) = [x(j1, j2), y(j1, j2), z(j1, j2)];
                count_points_on_edge = count_points_on_edge + 1;
                
            % Mark the points inside the area
            else
                points_in_element(count_points_in_element, :) = [x(j1, j2), y(j1, j2), z(j1, j2)];
                count_points_in_element = count_points_in_element + 1;
                
            end
        end
    end
    
    % Draw the points
    plot3(points_on_edge(:, 1), points_on_edge(:, 2), points_on_edge(:, 3), ...
          's', 'MarkerEdgeColor', [0.25; 0.25; 0.25], 'MarkerFaceColor', [0.20; 0.50; 0.80], 'MarkerSize', 8); hold on;
    
    plot3(points_in_element(:, 1), points_in_element(:, 2), points_in_element(:, 3), ...
          'h', 'MarkerEdgeColor', [0.30; 0.90; 0.25], 'MarkerFaceColor', [0.30; 0.90; 0.25], 'MarkerSize', 6); hold on;
    
    clear points_on_edge points_in_element;
    
    
    %----------------------------------------------------------------------
    %  Draw the surfaces
    %----------------------------------------------------------------------
    % Draw the left surface
    j1 = 1;
    plot3(x(j1, :), y(j1, :), z(j1, :), ...
          '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
    
    % Draw the right surface
    j1 = numPointsPerElementCoordinate;
    plot3(x(j1, :), y(j1, :), z(j1, :), ...
          '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
    
    % Draw the bottom surface
    j2 = 1;
    plot3(x(:, j2), y(:, j2), z(:, j2), ...
          '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
    
    % Draw the top surface
    j2 = numPointsPerElementCoordinate;
    plot3(x(:, j2), y(:, j2), z(:, j2), ...
          '-', 'LineWidth', 2.5, 'Color', [0.70; 0.35; 0.45]); hold on;
end