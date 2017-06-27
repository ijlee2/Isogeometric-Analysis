%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the 0th and 1st partial derivatives of 2D NURBS
%  with respect to the parametric domain at the given locations (e.g. at
%  the quadrature points).
%  
%  
%  Warning:
%  
%  We must first evaluate the 2D B-splines with respect to the parametric
%  domain. (Evaluate the univariate B-splines using build_bezier_extraction
%  and eval_1d_bernstein_der routines, and take the kronecker product.)
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [NURBS_der00, NURBS_der10, NURBS_der01] = eval_2d_nurbs_der1(w_e, Bspline_der00, Bspline_der10, Bspline_der01);
%  
%  where,
%  
%      k is the order up to which we differentiate NURBS
%      w_e is the weights for the NURBS nodes
%      Bspline_der00, Bspline_der10, Bspline_der01 are arrays of 2D
%          B-spline evaluations (size of numNodes x numQuadraturePoints)
%  
%  
%  Output:
%  
%  1. numNodes x numQuadraturePoints matrices
%  
%      The first column corresponds to the derivatives at the first point,
%      the second column to those at the second point, and so on.
%--------------------------------------------------------------------------
function [NURBS_der00, NURBS_der10, NURBS_der01] = eval_2d_nurbs_der1(w_e, Bspline_der00, Bspline_der10, Bspline_der01)
    % numNodes tells us how many basis functions there are on the element,
    % while numQuadraturePoints tells us how many evaluations we must make
    [numNodes, numQuadraturePoints] = size(Bspline_der00);
    
    
    % Initialize the arrays
    NURBS_der00 = zeros(numNodes, numQuadraturePoints);
    NURBS_der10 = NURBS_der00;
    NURBS_der01 = NURBS_der00;
    
    % Multiply the B-splines by the nodal weights
    for j = 1 : numQuadraturePoints
        Bspline_der00(:, j) = w_e .* Bspline_der00(:, j);
        Bspline_der10(:, j) = w_e .* Bspline_der10(:, j);
        Bspline_der01(:, j) = w_e .* Bspline_der01(:, j);
    end
    
    % Evaluate the 0th derivative
    sum00_inv = 1 ./ sum(Bspline_der00);
    for j = 1 : numQuadraturePoints
        NURBS_der00(:, j) = sum00_inv(j) * Bspline_der00(:, j);
    end
    
    % Evaluate the 1st partial derivatives
    sum10 = sum(Bspline_der10);
    sum01 = sum(Bspline_der01);
    for j = 1 : numQuadraturePoints
        NURBS_der10(:, j) = sum00_inv(j) * (Bspline_der10(:, j) - sum10(j) * NURBS_der00(:, j));
        NURBS_der01(:, j) = sum00_inv(j) * (Bspline_der01(:, j) - sum01(j) * NURBS_der00(:, j));
    end
end