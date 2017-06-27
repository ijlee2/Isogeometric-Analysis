%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the 0th and 1st partial derivatives of 3D NURBS
%  with respect to the parametric domain at the given locations (e.g. at
%  the quadrature points).
%  
%  
%  Warning:
%  
%  We must first evaluate the 3D B-splines with respect to the parametric
%  domain. (Evaluate the univariate B-splines using build_bezier_extraction
%  and eval_1d_bernstein_der routines, and take the kronecker product.)
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [NURBS_der000, NURBS_der100, NURBS_der010, NURBS_der001] = eval_3d_nurbs_der1(w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001);
%  
%  where,
%  
%      k is the order up to which we differentiate NURBS
%      w_e is the weights for the NURBS nodes
%      Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001 are
%          arrays of 3D B-spline evaluations (size of numNodes x
%          numQuadraturePoints)
%  
%  
%  Output:
%  
%  1. numNodes x numQuadraturePoints matrices
%  
%      The first column corresponds to the derivatives at the first point,
%      the second column to those at the second point, and so on.
%--------------------------------------------------------------------------
function [NURBS_der000, NURBS_der100, NURBS_der010, NURBS_der001] = eval_3d_nurbs_der1(w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001)
    % numNodes tells us how many basis functions there are on the element,
    % while numQuadraturePoints tells us how many evaluations we must make
    [numNodes, numQuadraturePoints] = size(Bspline_der000);
    
    
    % Initialize the arrays
    NURBS_der000 = zeros(numNodes, numQuadraturePoints);
    NURBS_der100 = NURBS_der000;
    NURBS_der010 = NURBS_der000;
    NURBS_der001 = NURBS_der000;
    
    % Multiply the B-splines by the nodal weights
    for j = 1 : numQuadraturePoints
        Bspline_der000(:, j) = w_e .* Bspline_der000(:, j);
        Bspline_der100(:, j) = w_e .* Bspline_der100(:, j);
        Bspline_der010(:, j) = w_e .* Bspline_der010(:, j);
        Bspline_der001(:, j) = w_e .* Bspline_der001(:, j);
    end
    
    % Evaluate the 0th derivative
    sum000_inv = 1 ./ sum(Bspline_der000);
    for j = 1 : numQuadraturePoints
        NURBS_der000(:, j) = sum000_inv(j) * Bspline_der000(:, j);
    end
    
    % Evaluate the 1st partial derivatives
    sum100 = sum(Bspline_der100);
    sum010 = sum(Bspline_der010);
    sum001 = sum(Bspline_der001);
    for j = 1 : numQuadraturePoints
        NURBS_der100(:, j) = sum000_inv(j) * (Bspline_der100(:, j) - sum100(j) * NURBS_der000(:, j));
        NURBS_der010(:, j) = sum000_inv(j) * (Bspline_der010(:, j) - sum010(j) * NURBS_der000(:, j));
        NURBS_der001(:, j) = sum000_inv(j) * (Bspline_der001(:, j) - sum001(j) * NURBS_der000(:, j));
    end
end