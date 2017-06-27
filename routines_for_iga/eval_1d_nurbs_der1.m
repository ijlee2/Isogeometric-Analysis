%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the 0th and 1st derivatives of 1D NURBS with
%  respect to the parametric domain at the given locations (e.g. at the
%  quadrature points).
%  
%  
%  Warning:
%  
%  We must first evaluate the 1D B-splines with respect to the parametric
%  domain. (Evaluate the 1D B-splines using build_bezier_extraction and
%  eval_1d_bernstein_der routines.)
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [NURBS_der0, NURBS_der1] = eval_1d_nurbs_der1(w_e, Bspline_der0, Bspline_der1);
%  
%  where,
%  
%      k is the order up to which we differentiate NURBS
%      w_e is the weights for the NURBS nodes
%      Bspline_der0, Bspline_der1 are arrays of 1D B-spline evaluations
%          (size of numNodes x numQuadraturePoints)
%  
%  
%  Output:
%  
%  1. numNodes x numQuadraturePoints matrices
%  
%      The first column corresponds to the derivatives at the first point,
%      the second column to those at the second point, and so on.
%--------------------------------------------------------------------------
function [NURBS_der0, NURBS_der1] = eval_1d_nurbs_der1(w_e, Bspline_der0, Bspline_der1)
    % numNodes tells us how many basis functions there are on the element,
    % while numQuadraturePoints tells us how many evaluations we must make
    [numNodes, numQuadraturePoints] = size(Bspline_der0);
    
    
    % Initialize the arrays
    NURBS_der0 = zeros(numNodes, numQuadraturePoints);
    NURBS_der1 = NURBS_der0;
    
    % Multiply the B-splines by the nodal weights
    for j = 1 : numQuadraturePoints
        Bspline_der0(:, j) = w_e .* Bspline_der0(:, j);
        Bspline_der1(:, j) = w_e .* Bspline_der1(:, j);
    end
    
    % Evaluate the 0th derivative
    sum0_inv = 1 ./ sum(Bspline_der0);
    for j = 1 : numQuadraturePoints
        NURBS_der0(:, j) = sum0_inv(j) * Bspline_der0(:, j);
    end
    
    % Evaluate the 1st derivative
    sum1 = sum(Bspline_der1);
    for j = 1 : numQuadraturePoints
        NURBS_der1(:, j) = sum0_inv(j) * (Bspline_der1(:, j) - sum1(j) * NURBS_der0(:, j));
    end
end