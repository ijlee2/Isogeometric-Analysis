%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine uses the de Casteljau algorithm to compute the Bernstein
%  polynomials of degree p at the given location.
%  
%  
%  Warning:
%  
%  The Bernstein domain is assumed to be [0, 1]. The coordinate is denoted
%  by t.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      output = eval_1d_bernstein(t, p);
%  
%  where,
%  
%      t is a point in the interval [0, 1]
%      p is the degree of the Bernstein polynomials
%  
%  
%  Output:
%  
%  1. (p + 1) x 1 column vector
%  
%      The i-th row returns the evaluation of the i-th Bernstein polynomial
%      at the location t.
%--------------------------------------------------------------------------
function output = eval_1d_bernstein(t, p)
    % The columns of B correspond to the coefficients for the Bernstein
    % polynomials basis functions
    B = eye(p + 1);
    
    % Use de Casteljau's algorithm
    for i = p : -1 : 1
        % Save memory by overwriting the entries of B
        for j = 1 : i
            B(:, j) = (1 - t) * B(:, j) + t * B(:, j + 1);
        end
    end
    
    output = B(:, 1);
end