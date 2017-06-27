%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine uses the de Casteljau algorithm to compute the derivatives
%  of the Bernstein polynomials of degree p at the given location.
%  
%  
%  Warning:
%  
%  The Bernstein domain is assumed to be [0, 1]. The coordinate is denoted
%  by t.
%  
%  If the derivatives of the Bernstein polynomials are not needed, please
%  use eval_bernstein_all instead.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      output = eval_1d_bernstein_der(t, p);
%  
%  where,
%  
%      t is a point in the interval [0, 1]
%      p is the degree of the Bernstein polynomials
%  
%  
%  Output:
%  
%  1. (p + 1) x (p + 1) matrix
%  
%      The first column corresponds to the 0th derivatives of the Bernstein
%      polynomials, the second column to the 1st derivatives, and so on,
%      at the location t.
%--------------------------------------------------------------------------
function output = eval_1d_bernstein_der(t, p)
    % Some useful constants
    constant_pp1 = p + 1;
    
    % The columns of B correspond to the coefficients for the Bernstein
    % polynomial basis functions
    B = [eye(constant_pp1); zeros(p, constant_pp1)];
    
    % Use de Casteljau's algorithm
    for i = p : -1 : 1
        % Save memory by overwriting the entries of B
        for j = 1 : i
            B(:, j) = (1 - t) * B(:, j) + t * B(:, j + 1);
        end
    end
    
    % Initialize the output matrix; the 0th derivative is already stored
    % in the matrix B
    output = zeros(constant_pp1);
    output(:, 1) = B(1 : constant_pp1, 1);
    
    % Compute the k-th derivative, up to the p-th derivative
    coefficient = 1;
    
    for k = 1 : p
        % Temporary variable
        constant_kp1 = k + 1;
        
        % Find the entries of B that are relevant in evaluting the k-th
        % derivative of the Bernstein polynomials
        temp = B(1 : (constant_pp1 + k), constant_kp1);
        
        % Compute the derivative using product rule
        sign = 1;
        
        for j = 1 : constant_kp1
            output(:, constant_kp1) = output(:, constant_kp1) + sign * nchoosek(k, j - 1) * temp(j : j + p);
            
            sign = -1 * sign;
        end
        
        % Multiply by the correct coefficient
        coefficient = coefficient * (p - k + 1);
        
        output(:, constant_kp1) = coefficient * output(:, constant_kp1);
    end
    
    %{
    %----------------------------------------------------------------------
    %  These are the matrices we expect to get for the quadratic and cubic
    %  Bernstein polynomials. Use them for debugging.
    %----------------------------------------------------------------------
    % True solution for the quadratic Bernstein polynomials
    truesol = [(1 - t)^2,       -2*(1 - t),         2;
               2*t*(1 - t),     2*(1 - 2*t),        -4;
               t^2,             2*t,                2];
    
    % True solution for the cubic Bernstein polynomials
    truesol = [(1 - t)^3,       -3*(1 - t)^2,           6*(1 - t),         -6;
               3*t*(1 - t)^2,   3*(1 - t)*(1 - 3*t),    6*(3*t - 2),       18;
               3*t^2*(1 - t),   -3*t*(3*t - 2),         6*(1 - 3*t),       -18;
               t^3,             3*t^2,                  6*t,               6];
    %}
end
