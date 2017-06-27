%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the derivatives of 2D NURBS with respect to the
%  parametric domain (up to the k-th partial derivative) at the given
%  locations (e.g. at the quadrature points).
%  
%  
%  Warning:
%  
%  We must first evaluate the 2D B-splines with respect to the parametric
%  domain. (Evaluate the univariate B-splines using build_bezier_extraction
%  and eval_1d_bernstein_der routines, and take the kronecker product.)
%  
%  For partial derivatives of order higher than 4, please see the general
%  formula in the official documentation. (A routine to generate the code
%  is included below.)
%  
%  If only the 0th and 1st derivatives are needed (e.g. for elasticity),
%  please use eval_2d_nurbs_der1 instead.
%  
%  
%  Instructions:
%  
%  Type one of the following onto Matlab's command window or in a code,
%  
%      [NURBS_der00, NURBS_der10, NURBS_der01, ... , NURBS_der04] = eval_2d_nurbs_derk(4, w_e, Bspline_der00, Bspline_der10, Bspline_der01, ... , Bspline_der04);
%      [NURBS_der00, NURBS_der10, NURBS_der01, ... , NURBS_der03] = eval_2d_nurbs_derk(3, w_e, Bspline_der00, Bspline_der10, Bspline_der01, ... , Bspline_der03);
%      [NURBS_der00, NURBS_der10, NURBS_der01, ... , NURBS_der02] = eval_2d_nurbs_derk(2, w_e, Bspline_der00, Bspline_der10, Bspline_der01, ... , Bspline_der02);
%  
%  where,
%  
%      k is the order up to which we differentiate NURBS
%      w_e is the weights for the NURBS nodes
%      Bspline_der00, etc. are arrays of 2D B-spline evaluations
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
function [NURBS_der00, ...
          NURBS_der10, NURBS_der01, ...
          NURBS_der20, NURBS_der11, NURBS_der02, ...
          NURBS_der30, NURBS_der21, NURBS_der12, NURBS_der03, ...
          NURBS_der40, NURBS_der31, NURBS_der22, NURBS_der13, NURBS_der04] = ...
eval_2d_nurbs_derk(k, w_e, ...
          Bspline_der00, ...
          Bspline_der10, Bspline_der01, ...
          Bspline_der20, Bspline_der11, Bspline_der02, ...
          Bspline_der30, Bspline_der21, Bspline_der12, Bspline_der03, ...
          Bspline_der40, Bspline_der31, Bspline_der22, Bspline_der13, Bspline_der04)
    % numNodes tells us how many basis functions there are on the element,
    % while numQuadraturePoints tells us how many evaluations we must make
    [numNodes, numQuadraturePoints] = size(Bspline_der00);
    
    
    % Compute the 0th derivative
    if (k == 0)
        % Initialize the arrays
        NURBS_der00 = zeros(numNodes, numQuadraturePoints);
        NURBS_der10 = [];
        NURBS_der01 = [];
        NURBS_der20 = [];
        NURBS_der11 = [];
        NURBS_der02 = [];
        NURBS_der30 = [];
        NURBS_der21 = [];
        NURBS_der12 = [];
        NURBS_der03 = [];
        NURBS_der40 = [];
        NURBS_der31 = [];
        NURBS_der22 = [];
        NURBS_der13 = [];
        NURBS_der04 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der00(:, j) = w_e .* Bspline_der00(:, j);
        end
        
        % Evaluate the 0th derivative
        sum00_inv = 1 ./ sum(Bspline_der00);
        for j = 1 : numQuadraturePoints
            NURBS_der00(:, j) = sum00_inv(j) * Bspline_der00(:, j);
        end
        
        
    % Compute the 0th and 1st derivatives
    elseif (k == 1)
        % Initialize the arrays
        NURBS_der00 = zeros(numNodes, numQuadraturePoints);
        NURBS_der10 = NURBS_der00;
        NURBS_der01 = NURBS_der00;
        NURBS_der20 = [];
        NURBS_der11 = [];
        NURBS_der02 = [];
        NURBS_der30 = [];
        NURBS_der21 = [];
        NURBS_der12 = [];
        NURBS_der03 = [];
        NURBS_der40 = [];
        NURBS_der31 = [];
        NURBS_der22 = [];
        NURBS_der13 = [];
        NURBS_der04 = [];
        
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
        
        
    % Compute the 0th, 1st, and 2nd derivatives
    elseif (k == 2)
        % Initialize the arrays
        NURBS_der00 = zeros(numNodes, numQuadraturePoints);
        NURBS_der10 = NURBS_der00;
        NURBS_der01 = NURBS_der00;
        NURBS_der20 = NURBS_der00;
        NURBS_der11 = NURBS_der00;
        NURBS_der02 = NURBS_der00;
        NURBS_der30 = [];
        NURBS_der21 = [];
        NURBS_der12 = [];
        NURBS_der03 = [];
        NURBS_der40 = [];
        NURBS_der31 = [];
        NURBS_der22 = [];
        NURBS_der13 = [];
        NURBS_der04 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der00(:, j) = w_e .* Bspline_der00(:, j);
            Bspline_der10(:, j) = w_e .* Bspline_der10(:, j);
            Bspline_der01(:, j) = w_e .* Bspline_der01(:, j);
            Bspline_der20(:, j) = w_e .* Bspline_der20(:, j);
            Bspline_der11(:, j) = w_e .* Bspline_der11(:, j);
            Bspline_der02(:, j) = w_e .* Bspline_der02(:, j);
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
        
        % Evaluate the 2nd partial derivatives
        sum20 = sum(Bspline_der20);
        sum11 = sum(Bspline_der11);
        sum02 = sum(Bspline_der02);
        for j = 1 : numQuadraturePoints
            NURBS_der20(:, j) = sum00_inv(j) * (Bspline_der20(:, j) - sum20(j) * NURBS_der00(:, j) - 2 * sum10(j) * NURBS_der10(:, j));
            NURBS_der11(:, j) = sum00_inv(j) * (Bspline_der11(:, j) - sum11(j) * NURBS_der00(:, j) -     sum01(j) * NURBS_der10(:, j) -     sum10(j) * NURBS_der01(:, j));
            NURBS_der02(:, j) = sum00_inv(j) * (Bspline_der02(:, j) - sum02(j) * NURBS_der00(:, j) - 2 * sum01(j) * NURBS_der01(:, j));
        end
        
        
    % Compute the 0th, 1st, 2nd, and 3rd derivatives
    elseif (k == 3)
        % Initialize the arrays
        NURBS_der00 = zeros(numNodes, numQuadraturePoints);
        NURBS_der10 = NURBS_der00;
        NURBS_der01 = NURBS_der00;
        NURBS_der20 = NURBS_der00;
        NURBS_der11 = NURBS_der00;
        NURBS_der02 = NURBS_der00;
        NURBS_der30 = NURBS_der00;
        NURBS_der21 = NURBS_der00;
        NURBS_der12 = NURBS_der00;
        NURBS_der03 = NURBS_der00;
        NURBS_der40 = [];
        NURBS_der31 = [];
        NURBS_der22 = [];
        NURBS_der13 = [];
        NURBS_der04 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der00(:, j) = w_e .* Bspline_der00(:, j);
            Bspline_der10(:, j) = w_e .* Bspline_der10(:, j);
            Bspline_der01(:, j) = w_e .* Bspline_der01(:, j);
            Bspline_der20(:, j) = w_e .* Bspline_der20(:, j);
            Bspline_der11(:, j) = w_e .* Bspline_der11(:, j);
            Bspline_der02(:, j) = w_e .* Bspline_der02(:, j);
            Bspline_der30(:, j) = w_e .* Bspline_der30(:, j);
            Bspline_der21(:, j) = w_e .* Bspline_der21(:, j);
            Bspline_der12(:, j) = w_e .* Bspline_der12(:, j);
            Bspline_der03(:, j) = w_e .* Bspline_der03(:, j);
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
        
        % Evaluate the 2nd partial derivatives
        sum20 = sum(Bspline_der20);
        sum11 = sum(Bspline_der11);
        sum02 = sum(Bspline_der02);
        for j = 1 : numQuadraturePoints
            NURBS_der20(:, j) = sum00_inv(j) * (Bspline_der20(:, j) - sum20(j) * NURBS_der00(:, j) - 2 * sum10(j) * NURBS_der10(:, j));
            NURBS_der11(:, j) = sum00_inv(j) * (Bspline_der11(:, j) - sum11(j) * NURBS_der00(:, j) -     sum01(j) * NURBS_der10(:, j) -     sum10(j) * NURBS_der01(:, j));
            NURBS_der02(:, j) = sum00_inv(j) * (Bspline_der02(:, j) - sum02(j) * NURBS_der00(:, j) - 2 * sum01(j) * NURBS_der01(:, j));
        end
        
        % Evaluate the 3rd partial derivatives
        sum30 = sum(Bspline_der30);
        sum21 = sum(Bspline_der21);
        sum12 = sum(Bspline_der12);
        sum03 = sum(Bspline_der03);
        for j = 1 : numQuadraturePoints
            NURBS_der30(:, j) = sum00_inv(j) * (Bspline_der30(:, j) - sum30(j) * NURBS_der00(:, j) - 3 * sum20(j) * NURBS_der10(:, j) - 3 * sum10(j) * NURBS_der20(:, j));
            NURBS_der21(:, j) = sum00_inv(j) * (Bspline_der21(:, j) - sum21(j) * NURBS_der00(:, j) - 2 * sum11(j) * NURBS_der10(:, j) -     sum01(j) * NURBS_der20(:, j) -     sum20(j) * NURBS_der01(:, j) - 2 * sum10(j) * NURBS_der11(:, j));
            NURBS_der12(:, j) = sum00_inv(j) * (Bspline_der12(:, j) - sum12(j) * NURBS_der00(:, j) -     sum02(j) * NURBS_der10(:, j) - 2 * sum11(j) * NURBS_der01(:, j) - 2 * sum01(j) * NURBS_der11(:, j) -     sum10(j) * NURBS_der02(:, j));
            NURBS_der03(:, j) = sum00_inv(j) * (Bspline_der03(:, j) - sum03(j) * NURBS_der00(:, j) - 3 * sum02(j) * NURBS_der01(:, j) - 3 * sum01(j) * NURBS_der02(:, j));
        end
        
        
    % Compute the 0th, 1st, 2nd, 3rd, and 4th derivatives
    elseif (k == 4)
        % Initialize the arrays
        NURBS_der00 = zeros(numNodes, numQuadraturePoints);
        NURBS_der10 = NURBS_der00;
        NURBS_der01 = NURBS_der00;
        NURBS_der20 = NURBS_der00;
        NURBS_der11 = NURBS_der00;
        NURBS_der02 = NURBS_der00;
        NURBS_der30 = NURBS_der00;
        NURBS_der21 = NURBS_der00;
        NURBS_der12 = NURBS_der00;
        NURBS_der03 = NURBS_der00;
        NURBS_der40 = NURBS_der00;
        NURBS_der31 = NURBS_der00;
        NURBS_der22 = NURBS_der00;
        NURBS_der13 = NURBS_der00;
        NURBS_der04 = NURBS_der00;
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der00(:, j) = w_e .* Bspline_der00(:, j);
            Bspline_der10(:, j) = w_e .* Bspline_der10(:, j);
            Bspline_der01(:, j) = w_e .* Bspline_der01(:, j);
            Bspline_der20(:, j) = w_e .* Bspline_der20(:, j);
            Bspline_der11(:, j) = w_e .* Bspline_der11(:, j);
            Bspline_der02(:, j) = w_e .* Bspline_der02(:, j);
            Bspline_der30(:, j) = w_e .* Bspline_der30(:, j);
            Bspline_der21(:, j) = w_e .* Bspline_der21(:, j);
            Bspline_der12(:, j) = w_e .* Bspline_der12(:, j);
            Bspline_der03(:, j) = w_e .* Bspline_der03(:, j);
            Bspline_der40(:, j) = w_e .* Bspline_der40(:, j);
            Bspline_der31(:, j) = w_e .* Bspline_der31(:, j);
            Bspline_der22(:, j) = w_e .* Bspline_der22(:, j);
            Bspline_der13(:, j) = w_e .* Bspline_der13(:, j);
            Bspline_der04(:, j) = w_e .* Bspline_der04(:, j);
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
        
        % Evaluate the 2nd partial derivatives
        sum20 = sum(Bspline_der20);
        sum11 = sum(Bspline_der11);
        sum02 = sum(Bspline_der02);
        for j = 1 : numQuadraturePoints
            NURBS_der20(:, j) = sum00_inv(j) * (Bspline_der20(:, j) - sum20(j) * NURBS_der00(:, j) - 2 * sum10(j) * NURBS_der10(:, j));
            NURBS_der11(:, j) = sum00_inv(j) * (Bspline_der11(:, j) - sum11(j) * NURBS_der00(:, j) -     sum01(j) * NURBS_der10(:, j) -     sum10(j) * NURBS_der01(:, j));
            NURBS_der02(:, j) = sum00_inv(j) * (Bspline_der02(:, j) - sum02(j) * NURBS_der00(:, j) - 2 * sum01(j) * NURBS_der01(:, j));
        end
        
        % Evaluate the 3rd partial derivatives
        sum30 = sum(Bspline_der30);
        sum21 = sum(Bspline_der21);
        sum12 = sum(Bspline_der12);
        sum03 = sum(Bspline_der03);
        for j = 1 : numQuadraturePoints
            NURBS_der30(:, j) = sum00_inv(j) * (Bspline_der30(:, j) - sum30(j) * NURBS_der00(:, j) - 3 * sum20(j) * NURBS_der10(:, j) - 3 * sum10(j) * NURBS_der20(:, j));
            NURBS_der21(:, j) = sum00_inv(j) * (Bspline_der21(:, j) - sum21(j) * NURBS_der00(:, j) - 2 * sum11(j) * NURBS_der10(:, j) -     sum01(j) * NURBS_der20(:, j) -     sum20(j) * NURBS_der01(:, j) - 2 * sum10(j) * NURBS_der11(:, j));
            NURBS_der12(:, j) = sum00_inv(j) * (Bspline_der12(:, j) - sum12(j) * NURBS_der00(:, j) -     sum02(j) * NURBS_der10(:, j) - 2 * sum11(j) * NURBS_der01(:, j) - 2 * sum01(j) * NURBS_der11(:, j) -     sum10(j) * NURBS_der02(:, j));
            NURBS_der03(:, j) = sum00_inv(j) * (Bspline_der03(:, j) - sum03(j) * NURBS_der00(:, j) - 3 * sum02(j) * NURBS_der01(:, j) - 3 * sum01(j) * NURBS_der02(:, j));
        end
        
        % Evaluate the 4th partial derivatives
        sum40 = sum(Bspline_der40);
        sum31 = sum(Bspline_der31);
        sum22 = sum(Bspline_der22);
        sum13 = sum(Bspline_der13);
        sum04 = sum(Bspline_der04);
        for j = 1 : numQuadraturePoints
            NURBS_der40(:, j) = sum00_inv(j) * (Bspline_der40(:, j) - sum40(j) * NURBS_der00(:, j) - 4 * sum30(j) * NURBS_der10(:, j) - 6 * sum20(j) * NURBS_der20(:, j) - 4 * sum10(j) * NURBS_der30(:, j));
            NURBS_der31(:, j) = sum00_inv(j) * (Bspline_der31(:, j) - sum31(j) * NURBS_der00(:, j) - 3 * sum21(j) * NURBS_der10(:, j) - 3 * sum11(j) * NURBS_der20(:, j) -     sum01(j) * NURBS_der30(:, j) -     sum30(j) * NURBS_der01(:, j) - 3 * sum20(j) * NURBS_der11(:, j) - 3 * sum10(j) * NURBS_der21(:, j));
            NURBS_der22(:, j) = sum00_inv(j) * (Bspline_der22(:, j) - sum22(j) * NURBS_der00(:, j) - 2 * sum12(j) * NURBS_der10(:, j) -     sum02(j) * NURBS_der20(:, j) - 2 * sum21(j) * NURBS_der01(:, j) - 4 * sum11(j) * NURBS_der11(:, j) - 2 * sum01(j) * NURBS_der21(:, j) -     sum20(j) * NURBS_der02(:, j) - 2 * sum10(j) * NURBS_der12(:, j));
            NURBS_der13(:, j) = sum00_inv(j) * (Bspline_der13(:, j) - sum13(j) * NURBS_der00(:, j) -     sum03(j) * NURBS_der10(:, j) - 3 * sum12(j) * NURBS_der01(:, j) - 3 * sum02(j) * NURBS_der11(:, j) - 3 * sum11(j) * NURBS_der02(:, j) - 3 * sum01(j) * NURBS_der12(:, j) -     sum10(j) * NURBS_der03(:, j));
            NURBS_der04(:, j) = sum00_inv(j) * (Bspline_der04(:, j) - sum04(j) * NURBS_der00(:, j) - 4 * sum03(j) * NURBS_der01(:, j) - 6 * sum02(j) * NURBS_der02(:, j) - 4 * sum01(j) * NURBS_der03(:, j));
        end
        
        
    end
end


%{
function code = code_to_generate_coefficients(k1, k2)
    k = k1 + k2;
    
    code = sprintf('NURBS_der%d%d(:, j) = sum00_inv(j) * (Bspline_der%d%d(:, j)', k1, k2, k1, k2);
    
    for j2 = 0 : min(k - 1, k2)
        for j1 = 0 : min(k - 1 - j2, k1)
            coefficient = nchoosek(k1, j1) * nchoosek(k2, j2);
            
            if (coefficient == 1)
                if (j1 == 0 && j2 == 0)
                    code = sprintf('%s - sum%d%d(j) * NURBS_der%d%d(:, j)', code, k1 - j1, k2 - j2, j1, j2);
                else
                    code = sprintf('%s -     sum%d%d(j) * NURBS_der%d%d(:, j)', code, k1 - j1, k2 - j2, j1, j2);
                end
            else
                code = sprintf('%s - %d * sum%d%d(j) * NURBS_der%d%d(:, j)', code, coefficient, k1 - j1, k2 - j2, j1, j2);
            end
        end
    end
    
    code = sprintf('%s);', code);
end
%}