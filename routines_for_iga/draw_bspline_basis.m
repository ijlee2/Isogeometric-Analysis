%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine draws the 1D B-splines of degree p or their derivatives
%  (up to the p-th derivative) for a given knot vector.
%  
%  
%  Warning:
%  
%  Here, the knot vector does not need to be open. The nodes are assumed
%  to be ordered along the direction of 1.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      draw_bspline_basis(knots, p, k);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      p is the degree of the B-splines
%      k is the order of the derivative
%  
%  
%  Output:
%  
%  1. Plot of the B-splines or their derivatives
%--------------------------------------------------------------------------
function draw_bspline_basis(knots, p, k)
    % Some useful constants
    constant_pp1 = p + 1;
    
    numKnots = size(knots, 1);
    
    % Check for errors
    if (k > p)
        fprintf('Error: The order of the derivative must be less than or equal to %d.\n\n', p);
        
        return;
    end
    
    
    %----------------------------------------------------------------------
    %  If the knot vector is not open, temporarily make it so
    %----------------------------------------------------------------------
    numKnotInsertions_begin = 0;
    for i = 1 : p
        if (knots(i) ~= knots(i + 1))
            numKnotInsertions_begin = constant_pp1 - i;
            
            break;
        end
    end
    
    numKnotInsertions_end = 0;
    for i = 1 : p
        if (knots(numKnots - i) ~= knots(numKnots - i + 1))
            numKnotInsertions_end = constant_pp1 - i;
            
            break;
        end
    end
    
    % Note that repmat correctly returns an empty column vector if the knot
    % vector is open at its end
    knots = [repmat(knots(1)       , numKnotInsertions_begin, 1); ...
             knots; ...
             repmat(knots(numKnots), numKnotInsertions_end  , 1)];
    numKnots = numKnots + numKnotInsertions_begin + numKnotInsertions_end;
    
    % Number of basis functions for the open knot vector
    numBasis = numKnots - constant_pp1;
    
    
    %----------------------------------------------------------------------
    %  Build the Bezier extraction matrices
    %----------------------------------------------------------------------
    [bezierExtractions, nodeIndexShifts, numElements] = build_bezier_extraction(knots, p);
    
    
    %----------------------------------------------------------------------
    %  Set parameters for plotting
    %----------------------------------------------------------------------
    figure;
    
    % Number of points where the fields are evaluated
    if (numElements > 8)
        numPointsPerElementCoordinate = 17;
    else
        numPointsPerElementCoordinate = 33;
    end
    
    % Points on the Bernstein domain
    t = linspace(0, 1, numPointsPerElementCoordinate)';
    
    % Evaluate the k-th derivatives of the Bernstein polynomials
    Bernstein_derk = zeros(constant_pp1, numPointsPerElementCoordinate);
    
    for j = 1 : numPointsPerElementCoordinate
        temp = eval_1d_bernstein_der(t(j), p);
        Bernstein_derk(:, j) = temp(:, k + 1);
    end
    
    % Initialize the points on the physical domain
    x = zeros(numPointsPerElementCoordinate, 1);
    y = zeros(numPointsPerElementCoordinate, numBasis);
    
    
    %----------------------------------------------------------------------
    %  Loop over the elements
    %----------------------------------------------------------------------
    for e = 1 : numElements
        % Get the positions of the knots
        xi1 = knots(e + p + nodeIndexShifts(e));
        xi2 = knots(e + constant_pp1 + nodeIndexShifts(e));
        
        % Find the Bezier extraction matrix
        bezierExtractions_e = bezierExtractions(:, :, e);
        
        
        %------------------------------------------------------------------
        %  Loop over the points
        %------------------------------------------------------------------
        % Evaluate dxi/dt at the points (constant)
        dxi_dt = xi2 - xi1;
        dt_dxi = 1 / dxi_dt;
        
        % Evaluate the k-th derivatives of the B-splines
        Bspline_derk = (dt_dxi^k * bezierExtractions_e) * Bernstein_derk;
        
        % Select knots uniformly from the B-spline domain [xi1, xi2]
        Delta_xi = dxi_dt / (numPointsPerElementCoordinate - 1);
        
        % Clear the entries
        y(:, :) = 0;
        
        for j = 1 : numPointsPerElementCoordinate
            % Find the corresponding point on the B-spline domain
            x(j) = xi1 + (j - 1) * Delta_xi;
            
            % Assign the basis function evaluations to the correct row
            y(j, (e : (e + p))' + nodeIndexShifts(e)) = Bspline_derk((1 : constant_pp1)', j);
        end
        
        
        %------------------------------------------------------------------
        %  Draw the B-splines
        %------------------------------------------------------------------
        for i = (1 + numKnotInsertions_begin) : (numBasis - numKnotInsertions_end)
            plot(x, y(:, i), '-', 'LineWidth', 2, 'Color', customColorMap((i - 1) / (numBasis - 1))); hold on;
        end
    end
    
    grid on;
    set(gca, 'FontSize', 18);
    if (k == 0)
        axis([knots(1) knots(numKnots) -0.05 1.05]);
        set(gca, 'YTick', linspace(0, 1, 11)');
    else
        axis([knots(1) knots(numKnots) -0.05 1.05]);
        axis 'auto y';
    end
end

% x = 0 --> red
% x = 1 --> cyan
function output = customColorMap(x)
    output = (1 - x) * [1; 0; 0] + x * [0; 1; 1];
end