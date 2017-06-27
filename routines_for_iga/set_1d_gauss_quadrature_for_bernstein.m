%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine returns the 1D Gauss quadrature points and weights for the
%  Bernstein domain [a, b] = [0, 1].
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [z, w] = set_1d_gauss_quadrature_for_bernstein(q);
%  
%  where,
%  
%      q is the number of quadrature points (degree of the quadrature)
%  
%  
%  Outputs:
%  
%  1. An array of Gauss points (in increasing nodal order)
%  
%  2. An array of Gauss weights
%--------------------------------------------------------------------------
function [z, w] = set_1d_gauss_quadrature_for_bernstein(q)
    % Find the quadrature rule for [a, b] = [0, 1]
    [z, w] = set_1d_gauss_quadrature(0, 1, q);
end