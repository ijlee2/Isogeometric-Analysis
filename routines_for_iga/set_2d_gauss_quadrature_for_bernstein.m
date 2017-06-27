%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine returns the 2D Gauss quadrature points and weights for the
%  Bernstein domain [0, 1] x [0, 1]. Unlike set_2d_gauss_quadrature, this
%  routine returns the arrays of 1D Gauss points, so that we can evaluate
%  the 2D B-splines efficiently as the kronecker product of 1D B-splines.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [z1, z2, w] = set_2d_gauss_quadrature_for_bernstein(q);
%  
%  where,
%  
%      q is an array of the number of quadrature points (column vector)
%  
%  
%  Outputs:
%  
%  1. Arrays of Gauss points (in increasing nodal order)
%  
%  2. An array of Gauss weights
%--------------------------------------------------------------------------
function [z1, z2, w] = set_2d_gauss_quadrature_for_bernstein(q)
    % Find the quadrature rule for [a1, b1] = [0, 1]
    [z1, w1] = set_1d_gauss_quadrature(0, 1, q(1));
    
    % Find the quadrature rule for [a2, b2] = [0, 1]
    [z2, w2] = set_1d_gauss_quadrature(0, 1, q(2));
    
    % Compute the weights for [a1, b1] x [a2, b2]
    w = kron(w2, w1);
end