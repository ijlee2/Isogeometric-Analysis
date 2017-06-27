%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the derivatives of 1D NURBS with respect to the
%  parametric domain (up to the k-th derivative) at the given locations
%  (e.g. at the quadrature points).
%  
%  
%  Warning:
%  
%  We must first evaluate the 1D B-splines with respect to the parametric
%  domain. (Evaluate the 1D B-splines using build_bezier_extraction and
%  eval_1d_bernstein_der routines.)
%  
%  For derivatives of order higher than 4, please see the general formula
%  in the official documentation. (A routine to generate the code is
%  included below.)
%  
%  If only the 0th and 1st derivatives are needed (e.g. for elasticity),
%  please use eval_1d_nurbs_der1 instead.
%  
%  
%  Instructions:
%  
%  Type one of the following onto Matlab's command window or in a code,
%  
%      [NURBS_der0, NURBS_der1, NURBS_der2, NURBS_der3, NURBS_der4] = eval_1d_nurbs_derk(4, w_e, Bspline_der0, Bspline_der1, Bspline_der2, Bspline_der3, Bspline_der4);
%      [NURBS_der0, NURBS_der1, NURBS_der2, NURBS_der3] = eval_1d_nurbs_derk(3, w_e, Bspline_der0, Bspline_der1, Bspline_der2, Bspline_der3);
%      [NURBS_der0, NURBS_der1, NURBS_der2] = eval_1d_nurbs_derk(2, w_e, Bspline_der0, Bspline_der1, Bspline_der2);
%  
%  where,
%  
%      k is the order up to which we differentiate NURBS
%      w_e is the weights for the NURBS nodes
%      Bspline_der0, etc. are arrays of 1D B-spline evaluations
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
function [NURBS_der0, NURBS_der1, NURBS_der2, NURBS_der3, NURBS_der4] = eval_1d_nurbs_derk(k, w_e, Bspline_der0, Bspline_der1, Bspline_der2, Bspline_der3, Bspline_der4)
    % numNodes tells us how many basis functions there are on the element,
    % while numQuadraturePoints tells us how many evaluations we must make
    [numNodes, numQuadraturePoints] = size(Bspline_der0);
    
    
    % Compute the 0th derivative
    if (k == 0)
        % Initialize the arrays
        NURBS_der0 = zeros(numNodes, numQuadraturePoints);
        NURBS_der1 = [];
        NURBS_der2 = [];
        NURBS_der3 = [];
        NURBS_der4 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der0(:, j) = w_e .* Bspline_der0(:, j);
        end
        
        % Evaluate the 0th derivative
        sum0_inv = 1 ./ sum(Bspline_der0);
        for j = 1 : numQuadraturePoints
            NURBS_der0(:, j) = sum0_inv(j) * Bspline_der0(:, j);
        end
        
        
    % Compute the 0th and 1st derivatives
    elseif (k == 1)
        % Initialize the arrays
        NURBS_der0 = zeros(numNodes, numQuadraturePoints);
        NURBS_der1 = NURBS_der0;
        NURBS_der2 = [];
        NURBS_der3 = [];
        NURBS_der4 = [];
        
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
        
        
    % Compute the 0th, 1st, and 2nd derivatives
    elseif (k == 2)
        % Initialize the arrays
        NURBS_der0 = zeros(numNodes, numQuadraturePoints);
        NURBS_der1 = NURBS_der0;
        NURBS_der2 = NURBS_der0;
        NURBS_der3 = [];
        NURBS_der4 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der0(:, j) = w_e .* Bspline_der0(:, j);
            Bspline_der1(:, j) = w_e .* Bspline_der1(:, j);
            Bspline_der2(:, j) = w_e .* Bspline_der2(:, j);
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
        
        % Evaluate the 2nd derivative
        sum2 = sum(Bspline_der2);
        for j = 1 : numQuadraturePoints
            NURBS_der2(:, j) = sum0_inv(j) * (Bspline_der2(:, j) - sum2(j) * NURBS_der0(:, j) - 2 * sum1(j) * NURBS_der1(:, j));
        end
        
        
    % Compute the 0th, 1st, 2nd, and 3rd derivatives
    elseif (k == 3)
        % Initialize the arrays
        NURBS_der0 = zeros(numNodes, numQuadraturePoints);
        NURBS_der1 = NURBS_der0;
        NURBS_der2 = NURBS_der0;
        NURBS_der3 = NURBS_der0;
        NURBS_der4 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der0(:, j) = w_e .* Bspline_der0(:, j);
            Bspline_der1(:, j) = w_e .* Bspline_der1(:, j);
            Bspline_der2(:, j) = w_e .* Bspline_der2(:, j);
            Bspline_der3(:, j) = w_e .* Bspline_der3(:, j);
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
        
        % Evaluate the 2nd derivative
        sum2 = sum(Bspline_der2);
        for j = 1 : numQuadraturePoints
            NURBS_der2(:, j) = sum0_inv(j) * (Bspline_der2(:, j) - sum2(j) * NURBS_der0(:, j) - 2 * sum1(j) * NURBS_der1(:, j));
        end
        
        % Evaluate the 3rd derivative
        sum3 = sum(Bspline_der3);
        for j = 1 : numQuadraturePoints
            NURBS_der3(:, j) = sum0_inv(j) * (Bspline_der3(:, j) - sum3(j) * NURBS_der0(:, j) - 3 * sum2(j) * NURBS_der1(:, j) - 3 * sum1(j) * NURBS_der2(:, j));
        end
        
        
    % Compute the 0th, 1st, 2nd, 3rd, and 4th derivatives
    elseif (k == 4)
        % Initialize the arrays
        NURBS_der0 = zeros(numNodes, numQuadraturePoints);
        NURBS_der1 = NURBS_der0;
        NURBS_der2 = NURBS_der0;
        NURBS_der3 = NURBS_der0;
        NURBS_der4 = NURBS_der0;
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der0(:, j) = w_e .* Bspline_der0(:, j);
            Bspline_der1(:, j) = w_e .* Bspline_der1(:, j);
            Bspline_der2(:, j) = w_e .* Bspline_der2(:, j);
            Bspline_der3(:, j) = w_e .* Bspline_der3(:, j);
            Bspline_der4(:, j) = w_e .* Bspline_der4(:, j);
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
        
        % Evaluate the 2nd derivative
        sum2 = sum(Bspline_der2);
        for j = 1 : numQuadraturePoints
            NURBS_der2(:, j) = sum0_inv(j) * (Bspline_der2(:, j) - sum2(j) * NURBS_der0(:, j) - 2 * sum1(j) * NURBS_der1(:, j));
        end
        
        % Evaluate the 3rd derivative
        sum3 = sum(Bspline_der3);
        for j = 1 : numQuadraturePoints
            NURBS_der3(:, j) = sum0_inv(j) * (Bspline_der3(:, j) - sum3(j) * NURBS_der0(:, j) - 3 * sum2(j) * NURBS_der1(:, j) - 3 * sum1(j) * NURBS_der2(:, j));
        end
        
        % Evaluate the 4th derivative
        sum4 = sum(Bspline_der4);
        for j = 1 : numQuadraturePoints
            NURBS_der3(:, j) = sum0_inv(j) * (Bspline_der4(:, j) - sum4(j) * NURBS_der0(:, j) - 4 * sum3(j) * NURBS_der1(:, j) - 6 * sum2(j) * NURBS_der2(:, j) - 4 * sum1(j) * NURBS_der3(:, j));
        end
        
        
    end
end


%{
function code = code_to_generate_coefficients(k)
    code = sprintf('NURBS_der%d(:, j) = sum0_inv(j) * (Bspline_der%d(:, j)', k, k);
    
    for j = 0 : k - 1
        coefficient = nchoosek(k, j);
        
        if (coefficient == 1)
            if (j == 0)
                code = sprintf('%s - sum%d(j) * NURBS_der%d(:, j)', code, k - j, j);
            else
                code = sprintf('%s -     sum%d(j) * NURBS_der%d(:, j)', code, k - j, j);
            end
        else
            code = sprintf('%s - %d * sum%d(j) * NURBS_der%d(:, j)', code, coefficient, k - j, j);
        end
    end
    
    code = sprintf('%s);', code);
end
%}