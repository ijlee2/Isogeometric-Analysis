%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine returns the 2D Gauss quadrature points and weights for the
%  domain [a1, b1] x [a2, b2].
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [z, w] = set_2d_gauss_quadrature(a, b, q);
%  
%  where,
%  
%      a is an array of the left endpoints (column vector)
%      b is an array of the right endpoints (column vector)
%      q is an array of the number of quadrature points (column vector)
%  
%  
%  Outputs:
%  
%  1. An array of Gauss points (in increasing nodal order)
%  
%  2. An array of Gauss weights
%--------------------------------------------------------------------------
function [z, w] = set_2d_gauss_quadrature(a, b, q)
    % Find the quadrature rule for [a1, b1]
    [z1, w1] = set_1d_gauss_quadrature(a(1), b(1), q(1));
    
    % Find the quadrature rule for [a2, b2]
    [z2, w2] = set_1d_gauss_quadrature(a(2), b(2), q(2));
    
    % Compute the quadrature points and weights for [a1, b1] x [a2, b2]
    z = [repmat(z1, q(2), 1), kron(z2, ones(q(1), 1))];
    w = kron(w2, w1);
end