%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine computes the derivatives of 3D NURBS with respect to the
%  parametric domain (up to the k-th partial derivative) at the given
%  locations (e.g. at the quadrature points).
%  
%  
%  Warning:
%  
%  We must first evaluate the 3D B-splines with respect to the parametric
%  domain. (Evaluate the univariate B-splines using build_bezier_extraction
%  and eval_1d_bernstein_der routines, and take the kronecker product.)
%  
%  For partial derivatives of order higher than 4, please see the general
%  formula in the official documentation. (A routine to generate the code
%  is included below.)
%  
%  If only the 0th and 1st derivatives are needed (e.g. for elasticity),
%  please use eval_3d_nurbs_der1 instead.
%  
%  
%  Instructions:
%  
%  Type one of the following onto Matlab's command window or in a code,
%  
%      [NURBS_der000, NURBS_der100, NURBS_der010, NURBS_der001, ... , NURBS_der004] = eval_3d_nurbs_derk(4, w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001, ... , Bspline_der004);
%      [NURBS_der000, NURBS_der100, NURBS_der010, NURBS_der001, ... , NURBS_der003] = eval_3d_nurbs_derk(3, w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001, ... , Bspline_der003);
%      [NURBS_der000, NURBS_der100, NURBS_der010, NURBS_der001, ... , NURBS_der002] = eval_3d_nurbs_derk(2, w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001, ... , Bspline_der002);
%  
%  where,
%  
%      k is the order up to which we differentiate NURBS
%      w_e is the weights for the NURBS nodes
%      Bspline_der000, etc. are arrays of 3D B-spline evaluations
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
function [NURBS_der000, ...
          NURBS_der100, NURBS_der010, NURBS_der001, ...
          NURBS_der200, NURBS_der110, NURBS_der101, NURBS_der020, NURBS_der011, NURBS_der002, ...
          NURBS_der300, NURBS_der210, NURBS_der201, NURBS_der120, NURBS_der111, NURBS_der102, NURBS_der030, NURBS_der021, NURBS_der012, NURBS_der003, ...
          NURBS_der400, NURBS_der310, NURBS_der301, NURBS_der220, NURBS_der211, NURBS_der202, NURBS_der130, NURBS_der121, NURBS_der112, NURBS_der103, NURBS_der040, NURBS_der031, NURBS_der022, NURBS_der013, NURBS_der004] = ...
eval_3d_nurbs_derk(k, w_e, ...
          Bspline_der000, ...
          Bspline_der100, Bspline_der010, Bspline_der001, ...
          Bspline_der200, Bspline_der110, Bspline_der101, Bspline_der020, Bspline_der011, Bspline_der002, ...
          Bspline_der300, Bspline_der210, Bspline_der201, Bspline_der120, Bspline_der111, Bspline_der102, Bspline_der030, Bspline_der021, Bspline_der012, Bspline_der003, ...
          Bspline_der400, Bspline_der310, Bspline_der301, Bspline_der220, Bspline_der211, Bspline_der202, Bspline_der130, Bspline_der121, Bspline_der112, Bspline_der103, Bspline_der040, Bspline_der031, Bspline_der022, Bspline_der013, Bspline_der004)
    % numNodes tells us how many basis functions there are on the element,
    % while numQuadraturePoints tells us how many evaluations we must make
    [numNodes, numQuadraturePoints] = size(Bspline_der000);
    
    
    % Compute the 0th derivative
    if (k == 0)
        % Initialize the arrays
        NURBS_der000 = zeros(numNodes, numQuadraturePoints);
        NURBS_der100 = [];
        NURBS_der010 = [];
        NURBS_der001 = [];
        NURBS_der200 = [];
        NURBS_der110 = [];
        NURBS_der101 = [];
        NURBS_der020 = [];
        NURBS_der011 = [];
        NURBS_der002 = [];
        NURBS_der300 = [];
        NURBS_der210 = [];
        NURBS_der201 = [];
        NURBS_der120 = [];
        NURBS_der111 = [];
        NURBS_der102 = [];
        NURBS_der030 = [];
        NURBS_der021 = [];
        NURBS_der012 = [];
        NURBS_der003 = [];
        NURBS_der400 = [];
        NURBS_der310 = [];
        NURBS_der301 = [];
        NURBS_der220 = [];
        NURBS_der211 = [];
        NURBS_der202 = [];
        NURBS_der130 = [];
        NURBS_der121 = [];
        NURBS_der112 = [];
        NURBS_der103 = [];
        NURBS_der040 = [];
        NURBS_der031 = [];
        NURBS_der022 = [];
        NURBS_der013 = [];
        NURBS_der004 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der000(:, j) = w_e .* Bspline_der000(:, j);
        end
        
        % Evaluate the 0th derivative
        sum000_inv = 1 ./ sum(Bspline_der000);
        for j = 1 : numQuadraturePoints
            NURBS_der000(:, j) = sum000_inv(j) * Bspline_der000(:, j);
        end
        
        
    % Compute the 0th and 1st derivatives
    elseif (k == 1)
        % Initialize the arrays
        NURBS_der000 = zeros(numNodes, numQuadraturePoints);
        NURBS_der100 = NURBS_der000;
        NURBS_der010 = NURBS_der000;
        NURBS_der001 = NURBS_der000;
        NURBS_der200 = [];
        NURBS_der110 = [];
        NURBS_der101 = [];
        NURBS_der020 = [];
        NURBS_der011 = [];
        NURBS_der002 = [];
        NURBS_der300 = [];
        NURBS_der210 = [];
        NURBS_der201 = [];
        NURBS_der120 = [];
        NURBS_der111 = [];
        NURBS_der102 = [];
        NURBS_der030 = [];
        NURBS_der021 = [];
        NURBS_der012 = [];
        NURBS_der003 = [];
        NURBS_der400 = [];
        NURBS_der310 = [];
        NURBS_der301 = [];
        NURBS_der220 = [];
        NURBS_der211 = [];
        NURBS_der202 = [];
        NURBS_der130 = [];
        NURBS_der121 = [];
        NURBS_der112 = [];
        NURBS_der103 = [];
        NURBS_der040 = [];
        NURBS_der031 = [];
        NURBS_der022 = [];
        NURBS_der013 = [];
        NURBS_der004 = [];
        
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
        
        
    % Compute the 0th, 1st, and 2nd derivatives
    elseif (k == 2)
        % Initialize the arrays
        NURBS_der000 = zeros(numNodes, numQuadraturePoints);
        NURBS_der100 = NURBS_der000;
        NURBS_der010 = NURBS_der000;
        NURBS_der001 = NURBS_der000;
        NURBS_der200 = NURBS_der000;
        NURBS_der110 = NURBS_der000;
        NURBS_der101 = NURBS_der000;
        NURBS_der020 = NURBS_der000;
        NURBS_der011 = NURBS_der000;
        NURBS_der002 = NURBS_der000;
        NURBS_der300 = [];
        NURBS_der210 = [];
        NURBS_der201 = [];
        NURBS_der120 = [];
        NURBS_der111 = [];
        NURBS_der102 = [];
        NURBS_der030 = [];
        NURBS_der021 = [];
        NURBS_der012 = [];
        NURBS_der003 = [];
        NURBS_der400 = [];
        NURBS_der310 = [];
        NURBS_der301 = [];
        NURBS_der220 = [];
        NURBS_der211 = [];
        NURBS_der202 = [];
        NURBS_der130 = [];
        NURBS_der121 = [];
        NURBS_der112 = [];
        NURBS_der103 = [];
        NURBS_der040 = [];
        NURBS_der031 = [];
        NURBS_der022 = [];
        NURBS_der013 = [];
        NURBS_der004 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der000(:, j) = w_e .* Bspline_der000(:, j);
            Bspline_der100(:, j) = w_e .* Bspline_der100(:, j);
            Bspline_der010(:, j) = w_e .* Bspline_der010(:, j);
            Bspline_der001(:, j) = w_e .* Bspline_der001(:, j);
            Bspline_der200(:, j) = w_e .* Bspline_der200(:, j);
            Bspline_der110(:, j) = w_e .* Bspline_der110(:, j);
            Bspline_der101(:, j) = w_e .* Bspline_der101(:, j);
            Bspline_der020(:, j) = w_e .* Bspline_der020(:, j);
            Bspline_der011(:, j) = w_e .* Bspline_der011(:, j);
            Bspline_der002(:, j) = w_e .* Bspline_der002(:, j);
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
        
        % Evaluate the 2nd partial derivatives
        sum200 = sum(Bspline_der200);
        sum110 = sum(Bspline_der110);
        sum101 = sum(Bspline_der101);
        sum020 = sum(Bspline_der020);
        sum011 = sum(Bspline_der011);
        sum002 = sum(Bspline_der002);
        for j = 1 : numQuadraturePoints
            NURBS_der200(:, j) = sum000_inv(j) * (Bspline_der200(:, j) - sum200(j) * NURBS_der000(:, j) - 2 * sum100(j) * NURBS_der100(:, j));
            NURBS_der110(:, j) = sum000_inv(j) * (Bspline_der110(:, j) - sum110(j) * NURBS_der000(:, j) -     sum010(j) * NURBS_der100(:, j) -     sum100(j) * NURBS_der010(:, j));
            NURBS_der101(:, j) = sum000_inv(j) * (Bspline_der101(:, j) - sum101(j) * NURBS_der000(:, j) -     sum001(j) * NURBS_der100(:, j) -     sum100(j) * NURBS_der001(:, j));
            NURBS_der020(:, j) = sum000_inv(j) * (Bspline_der020(:, j) - sum020(j) * NURBS_der000(:, j) - 2 * sum010(j) * NURBS_der010(:, j));
            NURBS_der011(:, j) = sum000_inv(j) * (Bspline_der011(:, j) - sum011(j) * NURBS_der000(:, j) -     sum001(j) * NURBS_der010(:, j) -     sum010(j) * NURBS_der001(:, j));
            NURBS_der002(:, j) = sum000_inv(j) * (Bspline_der002(:, j) - sum002(j) * NURBS_der000(:, j) - 2 * sum001(j) * NURBS_der001(:, j));
        end
        
        
    % Compute the 0th, 1st, 2nd, and 3rd derivatives
    elseif (k == 3)
        % Initialize the arrays
        NURBS_der000 = zeros(numNodes, numQuadraturePoints);
        NURBS_der100 = NURBS_der000;
        NURBS_der010 = NURBS_der000;
        NURBS_der001 = NURBS_der000;
        NURBS_der200 = NURBS_der000;
        NURBS_der110 = NURBS_der000;
        NURBS_der101 = NURBS_der000;
        NURBS_der020 = NURBS_der000;
        NURBS_der011 = NURBS_der000;
        NURBS_der002 = NURBS_der000;
        NURBS_der300 = NURBS_der000;
        NURBS_der210 = NURBS_der000;
        NURBS_der201 = NURBS_der000;
        NURBS_der120 = NURBS_der000;
        NURBS_der111 = NURBS_der000;
        NURBS_der102 = NURBS_der000;
        NURBS_der030 = NURBS_der000;
        NURBS_der021 = NURBS_der000;
        NURBS_der012 = NURBS_der000;
        NURBS_der003 = NURBS_der000;
        NURBS_der400 = [];
        NURBS_der310 = [];
        NURBS_der301 = [];
        NURBS_der220 = [];
        NURBS_der211 = [];
        NURBS_der202 = [];
        NURBS_der130 = [];
        NURBS_der121 = [];
        NURBS_der112 = [];
        NURBS_der103 = [];
        NURBS_der040 = [];
        NURBS_der031 = [];
        NURBS_der022 = [];
        NURBS_der013 = [];
        NURBS_der004 = [];
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der000(:, j) = w_e .* Bspline_der000(:, j);
            Bspline_der100(:, j) = w_e .* Bspline_der100(:, j);
            Bspline_der010(:, j) = w_e .* Bspline_der010(:, j);
            Bspline_der001(:, j) = w_e .* Bspline_der001(:, j);
            Bspline_der200(:, j) = w_e .* Bspline_der200(:, j);
            Bspline_der110(:, j) = w_e .* Bspline_der110(:, j);
            Bspline_der101(:, j) = w_e .* Bspline_der101(:, j);
            Bspline_der020(:, j) = w_e .* Bspline_der020(:, j);
            Bspline_der011(:, j) = w_e .* Bspline_der011(:, j);
            Bspline_der002(:, j) = w_e .* Bspline_der002(:, j);
            Bspline_der300(:, j) = w_e .* Bspline_der300(:, j);
            Bspline_der210(:, j) = w_e .* Bspline_der210(:, j);
            Bspline_der201(:, j) = w_e .* Bspline_der201(:, j);
            Bspline_der120(:, j) = w_e .* Bspline_der120(:, j);
            Bspline_der111(:, j) = w_e .* Bspline_der111(:, j);
            Bspline_der102(:, j) = w_e .* Bspline_der102(:, j);
            Bspline_der030(:, j) = w_e .* Bspline_der030(:, j);
            Bspline_der021(:, j) = w_e .* Bspline_der021(:, j);
            Bspline_der012(:, j) = w_e .* Bspline_der012(:, j);
            Bspline_der003(:, j) = w_e .* Bspline_der003(:, j);
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
        
        % Evaluate the 2nd partial derivatives
        sum200 = sum(Bspline_der200);
        sum110 = sum(Bspline_der110);
        sum101 = sum(Bspline_der101);
        sum020 = sum(Bspline_der020);
        sum011 = sum(Bspline_der011);
        sum002 = sum(Bspline_der002);
        for j = 1 : numQuadraturePoints
            NURBS_der200(:, j) = sum000_inv(j) * (Bspline_der200(:, j) - sum200(j) * NURBS_der000(:, j) - 2 * sum100(j) * NURBS_der100(:, j));
            NURBS_der110(:, j) = sum000_inv(j) * (Bspline_der110(:, j) - sum110(j) * NURBS_der000(:, j) -     sum010(j) * NURBS_der100(:, j) -     sum100(j) * NURBS_der010(:, j));
            NURBS_der101(:, j) = sum000_inv(j) * (Bspline_der101(:, j) - sum101(j) * NURBS_der000(:, j) -     sum001(j) * NURBS_der100(:, j) -     sum100(j) * NURBS_der001(:, j));
            NURBS_der020(:, j) = sum000_inv(j) * (Bspline_der020(:, j) - sum020(j) * NURBS_der000(:, j) - 2 * sum010(j) * NURBS_der010(:, j));
            NURBS_der011(:, j) = sum000_inv(j) * (Bspline_der011(:, j) - sum011(j) * NURBS_der000(:, j) -     sum001(j) * NURBS_der010(:, j) -     sum010(j) * NURBS_der001(:, j));
            NURBS_der002(:, j) = sum000_inv(j) * (Bspline_der002(:, j) - sum002(j) * NURBS_der000(:, j) - 2 * sum001(j) * NURBS_der001(:, j));
        end
        
        % Evaluate the 3rd partial derivatives
        sum300 = sum(Bspline_der300);
        sum210 = sum(Bspline_der210);
        sum201 = sum(Bspline_der201);
        sum120 = sum(Bspline_der120);
        sum111 = sum(Bspline_der111);
        sum102 = sum(Bspline_der102);
        sum030 = sum(Bspline_der030);
        sum021 = sum(Bspline_der021);
        sum012 = sum(Bspline_der012);
        sum003 = sum(Bspline_der003);
        for j = 1 : numQuadraturePoints
            NURBS_der300(:, j) = sum000_inv(j) * (Bspline_der300(:, j) - sum300(j) * NURBS_der000(:, j) - 3 * sum200(j) * NURBS_der100(:, j) - 3 * sum100(j) * NURBS_der200(:, j));
            NURBS_der210(:, j) = sum000_inv(j) * (Bspline_der210(:, j) - sum210(j) * NURBS_der000(:, j) - 2 * sum110(j) * NURBS_der100(:, j) -     sum010(j) * NURBS_der200(:, j) -     sum200(j) * NURBS_der010(:, j) - 2 * sum100(j) * NURBS_der110(:, j));
            NURBS_der201(:, j) = sum000_inv(j) * (Bspline_der201(:, j) - sum201(j) * NURBS_der000(:, j) - 2 * sum101(j) * NURBS_der100(:, j) -     sum001(j) * NURBS_der200(:, j) -     sum200(j) * NURBS_der001(:, j) - 2 * sum100(j) * NURBS_der101(:, j));
            NURBS_der120(:, j) = sum000_inv(j) * (Bspline_der120(:, j) - sum120(j) * NURBS_der000(:, j) -     sum020(j) * NURBS_der100(:, j) - 2 * sum110(j) * NURBS_der010(:, j) - 2 * sum010(j) * NURBS_der110(:, j) -     sum100(j) * NURBS_der020(:, j));
            NURBS_der111(:, j) = sum000_inv(j) * (Bspline_der111(:, j) - sum111(j) * NURBS_der000(:, j) -     sum011(j) * NURBS_der100(:, j) -     sum101(j) * NURBS_der010(:, j) -     sum001(j) * NURBS_der110(:, j) -     sum110(j) * NURBS_der001(:, j) -     sum010(j) * NURBS_der101(:, j) -     sum100(j) * NURBS_der011(:, j));
            NURBS_der102(:, j) = sum000_inv(j) * (Bspline_der102(:, j) - sum102(j) * NURBS_der000(:, j) -     sum002(j) * NURBS_der100(:, j) - 2 * sum101(j) * NURBS_der001(:, j) - 2 * sum001(j) * NURBS_der101(:, j) -     sum100(j) * NURBS_der002(:, j));
            NURBS_der030(:, j) = sum000_inv(j) * (Bspline_der030(:, j) - sum030(j) * NURBS_der000(:, j) - 3 * sum020(j) * NURBS_der010(:, j) - 3 * sum010(j) * NURBS_der020(:, j));
            NURBS_der021(:, j) = sum000_inv(j) * (Bspline_der021(:, j) - sum021(j) * NURBS_der000(:, j) - 2 * sum011(j) * NURBS_der010(:, j) -     sum001(j) * NURBS_der020(:, j) -     sum020(j) * NURBS_der001(:, j) - 2 * sum010(j) * NURBS_der011(:, j));
            NURBS_der012(:, j) = sum000_inv(j) * (Bspline_der012(:, j) - sum012(j) * NURBS_der000(:, j) -     sum002(j) * NURBS_der010(:, j) - 2 * sum011(j) * NURBS_der001(:, j) - 2 * sum001(j) * NURBS_der011(:, j) -     sum010(j) * NURBS_der002(:, j));
            NURBS_der003(:, j) = sum000_inv(j) * (Bspline_der003(:, j) - sum003(j) * NURBS_der000(:, j) - 3 * sum002(j) * NURBS_der001(:, j) - 3 * sum001(j) * NURBS_der002(:, j));
        end
        
        
    % Compute the 0th, 1st, 2nd, 3rd, and 4th derivatives
    elseif (k == 4)
        % Initialize the arrays
        NURBS_der000 = zeros(numNodes, numQuadraturePoints);
        NURBS_der100 = NURBS_der000;
        NURBS_der010 = NURBS_der000;
        NURBS_der001 = NURBS_der000;
        NURBS_der200 = NURBS_der000;
        NURBS_der110 = NURBS_der000;
        NURBS_der101 = NURBS_der000;
        NURBS_der020 = NURBS_der000;
        NURBS_der011 = NURBS_der000;
        NURBS_der002 = NURBS_der000;
        NURBS_der300 = NURBS_der000;
        NURBS_der210 = NURBS_der000;
        NURBS_der201 = NURBS_der000;
        NURBS_der120 = NURBS_der000;
        NURBS_der111 = NURBS_der000;
        NURBS_der102 = NURBS_der000;
        NURBS_der030 = NURBS_der000;
        NURBS_der021 = NURBS_der000;
        NURBS_der012 = NURBS_der000;
        NURBS_der003 = NURBS_der000;
        NURBS_der400 = NURBS_der000;
        NURBS_der310 = NURBS_der000;
        NURBS_der301 = NURBS_der000;
        NURBS_der220 = NURBS_der000;
        NURBS_der211 = NURBS_der000;
        NURBS_der202 = NURBS_der000;
        NURBS_der130 = NURBS_der000;
        NURBS_der121 = NURBS_der000;
        NURBS_der112 = NURBS_der000;
        NURBS_der103 = NURBS_der000;
        NURBS_der040 = NURBS_der000;
        NURBS_der031 = NURBS_der000;
        NURBS_der022 = NURBS_der000;
        NURBS_der013 = NURBS_der000;
        NURBS_der004 = NURBS_der000;
        
        % Multiply the B-splines by the nodal weights
        for j = 1 : numQuadraturePoints
            Bspline_der000(:, j) = w_e .* Bspline_der000(:, j);
            Bspline_der100(:, j) = w_e .* Bspline_der100(:, j);
            Bspline_der010(:, j) = w_e .* Bspline_der010(:, j);
            Bspline_der001(:, j) = w_e .* Bspline_der001(:, j);
            Bspline_der200(:, j) = w_e .* Bspline_der200(:, j);
            Bspline_der110(:, j) = w_e .* Bspline_der110(:, j);
            Bspline_der101(:, j) = w_e .* Bspline_der101(:, j);
            Bspline_der020(:, j) = w_e .* Bspline_der020(:, j);
            Bspline_der011(:, j) = w_e .* Bspline_der011(:, j);
            Bspline_der002(:, j) = w_e .* Bspline_der002(:, j);
            Bspline_der300(:, j) = w_e .* Bspline_der300(:, j);
            Bspline_der210(:, j) = w_e .* Bspline_der210(:, j);
            Bspline_der201(:, j) = w_e .* Bspline_der201(:, j);
            Bspline_der120(:, j) = w_e .* Bspline_der120(:, j);
            Bspline_der111(:, j) = w_e .* Bspline_der111(:, j);
            Bspline_der102(:, j) = w_e .* Bspline_der102(:, j);
            Bspline_der030(:, j) = w_e .* Bspline_der030(:, j);
            Bspline_der021(:, j) = w_e .* Bspline_der021(:, j);
            Bspline_der012(:, j) = w_e .* Bspline_der012(:, j);
            Bspline_der003(:, j) = w_e .* Bspline_der003(:, j);
            Bspline_der400(:, j) = w_e .* Bspline_der400(:, j);
            Bspline_der310(:, j) = w_e .* Bspline_der310(:, j);
            Bspline_der301(:, j) = w_e .* Bspline_der301(:, j);
            Bspline_der220(:, j) = w_e .* Bspline_der220(:, j);
            Bspline_der211(:, j) = w_e .* Bspline_der211(:, j);
            Bspline_der202(:, j) = w_e .* Bspline_der202(:, j);
            Bspline_der130(:, j) = w_e .* Bspline_der130(:, j);
            Bspline_der121(:, j) = w_e .* Bspline_der121(:, j);
            Bspline_der112(:, j) = w_e .* Bspline_der112(:, j);
            Bspline_der103(:, j) = w_e .* Bspline_der103(:, j);
            Bspline_der040(:, j) = w_e .* Bspline_der040(:, j);
            Bspline_der031(:, j) = w_e .* Bspline_der031(:, j);
            Bspline_der022(:, j) = w_e .* Bspline_der022(:, j);
            Bspline_der013(:, j) = w_e .* Bspline_der013(:, j);
            Bspline_der004(:, j) = w_e .* Bspline_der004(:, j);
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
        
        % Evaluate the 2nd partial derivatives
        sum200 = sum(Bspline_der200);
        sum110 = sum(Bspline_der110);
        sum101 = sum(Bspline_der101);
        sum020 = sum(Bspline_der020);
        sum011 = sum(Bspline_der011);
        sum002 = sum(Bspline_der002);
        for j = 1 : numQuadraturePoints
            NURBS_der200(:, j) = sum000_inv(j) * (Bspline_der200(:, j) - sum200(j) * NURBS_der000(:, j) - 2 * sum100(j) * NURBS_der100(:, j));
            NURBS_der110(:, j) = sum000_inv(j) * (Bspline_der110(:, j) - sum110(j) * NURBS_der000(:, j) -     sum010(j) * NURBS_der100(:, j) -     sum100(j) * NURBS_der010(:, j));
            NURBS_der101(:, j) = sum000_inv(j) * (Bspline_der101(:, j) - sum101(j) * NURBS_der000(:, j) -     sum001(j) * NURBS_der100(:, j) -     sum100(j) * NURBS_der001(:, j));
            NURBS_der020(:, j) = sum000_inv(j) * (Bspline_der020(:, j) - sum020(j) * NURBS_der000(:, j) - 2 * sum010(j) * NURBS_der010(:, j));
            NURBS_der011(:, j) = sum000_inv(j) * (Bspline_der011(:, j) - sum011(j) * NURBS_der000(:, j) -     sum001(j) * NURBS_der010(:, j) -     sum010(j) * NURBS_der001(:, j));
            NURBS_der002(:, j) = sum000_inv(j) * (Bspline_der002(:, j) - sum002(j) * NURBS_der000(:, j) - 2 * sum001(j) * NURBS_der001(:, j));
        end
        
        % Evaluate the 3rd partial derivatives
        sum300 = sum(Bspline_der300);
        sum210 = sum(Bspline_der210);
        sum201 = sum(Bspline_der201);
        sum120 = sum(Bspline_der120);
        sum111 = sum(Bspline_der111);
        sum102 = sum(Bspline_der102);
        sum030 = sum(Bspline_der030);
        sum021 = sum(Bspline_der021);
        sum012 = sum(Bspline_der012);
        sum003 = sum(Bspline_der003);
        for j = 1 : numQuadraturePoints
            NURBS_der300(:, j) = sum000_inv(j) * (Bspline_der300(:, j) - sum300(j) * NURBS_der000(:, j) - 3 * sum200(j) * NURBS_der100(:, j) - 3 * sum100(j) * NURBS_der200(:, j));
            NURBS_der210(:, j) = sum000_inv(j) * (Bspline_der210(:, j) - sum210(j) * NURBS_der000(:, j) - 2 * sum110(j) * NURBS_der100(:, j) -     sum010(j) * NURBS_der200(:, j) -     sum200(j) * NURBS_der010(:, j) - 2 * sum100(j) * NURBS_der110(:, j));
            NURBS_der201(:, j) = sum000_inv(j) * (Bspline_der201(:, j) - sum201(j) * NURBS_der000(:, j) - 2 * sum101(j) * NURBS_der100(:, j) -     sum001(j) * NURBS_der200(:, j) -     sum200(j) * NURBS_der001(:, j) - 2 * sum100(j) * NURBS_der101(:, j));
            NURBS_der120(:, j) = sum000_inv(j) * (Bspline_der120(:, j) - sum120(j) * NURBS_der000(:, j) -     sum020(j) * NURBS_der100(:, j) - 2 * sum110(j) * NURBS_der010(:, j) - 2 * sum010(j) * NURBS_der110(:, j) -     sum100(j) * NURBS_der020(:, j));
            NURBS_der111(:, j) = sum000_inv(j) * (Bspline_der111(:, j) - sum111(j) * NURBS_der000(:, j) -     sum011(j) * NURBS_der100(:, j) -     sum101(j) * NURBS_der010(:, j) -     sum001(j) * NURBS_der110(:, j) -     sum110(j) * NURBS_der001(:, j) -     sum010(j) * NURBS_der101(:, j) -     sum100(j) * NURBS_der011(:, j));
            NURBS_der102(:, j) = sum000_inv(j) * (Bspline_der102(:, j) - sum102(j) * NURBS_der000(:, j) -     sum002(j) * NURBS_der100(:, j) - 2 * sum101(j) * NURBS_der001(:, j) - 2 * sum001(j) * NURBS_der101(:, j) -     sum100(j) * NURBS_der002(:, j));
            NURBS_der030(:, j) = sum000_inv(j) * (Bspline_der030(:, j) - sum030(j) * NURBS_der000(:, j) - 3 * sum020(j) * NURBS_der010(:, j) - 3 * sum010(j) * NURBS_der020(:, j));
            NURBS_der021(:, j) = sum000_inv(j) * (Bspline_der021(:, j) - sum021(j) * NURBS_der000(:, j) - 2 * sum011(j) * NURBS_der010(:, j) -     sum001(j) * NURBS_der020(:, j) -     sum020(j) * NURBS_der001(:, j) - 2 * sum010(j) * NURBS_der011(:, j));
            NURBS_der012(:, j) = sum000_inv(j) * (Bspline_der012(:, j) - sum012(j) * NURBS_der000(:, j) -     sum002(j) * NURBS_der010(:, j) - 2 * sum011(j) * NURBS_der001(:, j) - 2 * sum001(j) * NURBS_der011(:, j) -     sum010(j) * NURBS_der002(:, j));
            NURBS_der003(:, j) = sum000_inv(j) * (Bspline_der003(:, j) - sum003(j) * NURBS_der000(:, j) - 3 * sum002(j) * NURBS_der001(:, j) - 3 * sum001(j) * NURBS_der002(:, j));
        end
        
        % Evaluate the 4th partial derivatives
        sum400 = sum(Bspline_der400);
        sum310 = sum(Bspline_der310);
        sum301 = sum(Bspline_der301);
        sum220 = sum(Bspline_der220);
        sum211 = sum(Bspline_der211);
        sum202 = sum(Bspline_der202);
        sum130 = sum(Bspline_der130);
        sum121 = sum(Bspline_der121);
        sum112 = sum(Bspline_der112);
        sum103 = sum(Bspline_der103);
        sum040 = sum(Bspline_der040);
        sum031 = sum(Bspline_der031);
        sum022 = sum(Bspline_der022);
        sum013 = sum(Bspline_der013);
        sum004 = sum(Bspline_der004);
        for j = 1 : numQuadraturePoints
            NURBS_der400(:, j) = sum000_inv(j) * (Bspline_der400(:, j) - sum400(j) * NURBS_der000(:, j) - 4 * sum300(j) * NURBS_der100(:, j) - 6 * sum200(j) * NURBS_der200(:, j) - 4 * sum100(j) * NURBS_der300(:, j));
            NURBS_der310(:, j) = sum000_inv(j) * (Bspline_der310(:, j) - sum310(j) * NURBS_der000(:, j) - 3 * sum210(j) * NURBS_der100(:, j) - 3 * sum110(j) * NURBS_der200(:, j) -     sum010(j) * NURBS_der300(:, j) -     sum300(j) * NURBS_der010(:, j) - 3 * sum200(j) * NURBS_der110(:, j) - 3 * sum100(j) * NURBS_der210(:, j));
            NURBS_der301(:, j) = sum000_inv(j) * (Bspline_der301(:, j) - sum301(j) * NURBS_der000(:, j) - 3 * sum201(j) * NURBS_der100(:, j) - 3 * sum101(j) * NURBS_der200(:, j) -     sum001(j) * NURBS_der300(:, j) -     sum300(j) * NURBS_der001(:, j) - 3 * sum200(j) * NURBS_der101(:, j) - 3 * sum100(j) * NURBS_der201(:, j));
            NURBS_der220(:, j) = sum000_inv(j) * (Bspline_der220(:, j) - sum220(j) * NURBS_der000(:, j) - 2 * sum120(j) * NURBS_der100(:, j) -     sum020(j) * NURBS_der200(:, j) - 2 * sum210(j) * NURBS_der010(:, j) - 4 * sum110(j) * NURBS_der110(:, j) - 2 * sum010(j) * NURBS_der210(:, j) -     sum200(j) * NURBS_der020(:, j) - 2 * sum100(j) * NURBS_der120(:, j));
            NURBS_der211(:, j) = sum000_inv(j) * (Bspline_der211(:, j) - sum211(j) * NURBS_der000(:, j) - 2 * sum111(j) * NURBS_der100(:, j) -     sum011(j) * NURBS_der200(:, j) -     sum201(j) * NURBS_der010(:, j) - 2 * sum101(j) * NURBS_der110(:, j) -     sum001(j) * NURBS_der210(:, j) -     sum210(j) * NURBS_der001(:, j) - 2 * sum110(j) * NURBS_der101(:, j) -     sum010(j) * NURBS_der201(:, j) -     sum200(j) * NURBS_der011(:, j) - 2 * sum100(j) * NURBS_der111(:, j));
            NURBS_der202(:, j) = sum000_inv(j) * (Bspline_der202(:, j) - sum202(j) * NURBS_der000(:, j) - 2 * sum102(j) * NURBS_der100(:, j) -     sum002(j) * NURBS_der200(:, j) - 2 * sum201(j) * NURBS_der001(:, j) - 4 * sum101(j) * NURBS_der101(:, j) - 2 * sum001(j) * NURBS_der201(:, j) -     sum200(j) * NURBS_der002(:, j) - 2 * sum100(j) * NURBS_der102(:, j));
            NURBS_der130(:, j) = sum000_inv(j) * (Bspline_der130(:, j) - sum130(j) * NURBS_der000(:, j) -     sum030(j) * NURBS_der100(:, j) - 3 * sum120(j) * NURBS_der010(:, j) - 3 * sum020(j) * NURBS_der110(:, j) - 3 * sum110(j) * NURBS_der020(:, j) - 3 * sum010(j) * NURBS_der120(:, j) -     sum100(j) * NURBS_der030(:, j));
            NURBS_der121(:, j) = sum000_inv(j) * (Bspline_der121(:, j) - sum121(j) * NURBS_der000(:, j) -     sum021(j) * NURBS_der100(:, j) - 2 * sum111(j) * NURBS_der010(:, j) - 2 * sum011(j) * NURBS_der110(:, j) -     sum101(j) * NURBS_der020(:, j) -     sum001(j) * NURBS_der120(:, j) -     sum120(j) * NURBS_der001(:, j) -     sum020(j) * NURBS_der101(:, j) - 2 * sum110(j) * NURBS_der011(:, j) - 2 * sum010(j) * NURBS_der111(:, j) -     sum100(j) * NURBS_der021(:, j));
            NURBS_der112(:, j) = sum000_inv(j) * (Bspline_der112(:, j) - sum112(j) * NURBS_der000(:, j) -     sum012(j) * NURBS_der100(:, j) -     sum102(j) * NURBS_der010(:, j) -     sum002(j) * NURBS_der110(:, j) - 2 * sum111(j) * NURBS_der001(:, j) - 2 * sum011(j) * NURBS_der101(:, j) - 2 * sum101(j) * NURBS_der011(:, j) - 2 * sum001(j) * NURBS_der111(:, j) -     sum110(j) * NURBS_der002(:, j) -     sum010(j) * NURBS_der102(:, j) -     sum100(j) * NURBS_der012(:, j));
            NURBS_der103(:, j) = sum000_inv(j) * (Bspline_der103(:, j) - sum103(j) * NURBS_der000(:, j) -     sum003(j) * NURBS_der100(:, j) - 3 * sum102(j) * NURBS_der001(:, j) - 3 * sum002(j) * NURBS_der101(:, j) - 3 * sum101(j) * NURBS_der002(:, j) - 3 * sum001(j) * NURBS_der102(:, j) -     sum100(j) * NURBS_der003(:, j));
            NURBS_der040(:, j) = sum000_inv(j) * (Bspline_der040(:, j) - sum040(j) * NURBS_der000(:, j) - 4 * sum030(j) * NURBS_der010(:, j) - 6 * sum020(j) * NURBS_der020(:, j) - 4 * sum010(j) * NURBS_der030(:, j));
            NURBS_der031(:, j) = sum000_inv(j) * (Bspline_der031(:, j) - sum031(j) * NURBS_der000(:, j) - 3 * sum021(j) * NURBS_der010(:, j) - 3 * sum011(j) * NURBS_der020(:, j) -     sum001(j) * NURBS_der030(:, j) -     sum030(j) * NURBS_der001(:, j) - 3 * sum020(j) * NURBS_der011(:, j) - 3 * sum010(j) * NURBS_der021(:, j));
            NURBS_der022(:, j) = sum000_inv(j) * (Bspline_der022(:, j) - sum022(j) * NURBS_der000(:, j) - 2 * sum012(j) * NURBS_der010(:, j) -     sum002(j) * NURBS_der020(:, j) - 2 * sum021(j) * NURBS_der001(:, j) - 4 * sum011(j) * NURBS_der011(:, j) - 2 * sum001(j) * NURBS_der021(:, j) -     sum020(j) * NURBS_der002(:, j) - 2 * sum010(j) * NURBS_der012(:, j));
            NURBS_der013(:, j) = sum000_inv(j) * (Bspline_der013(:, j) - sum013(j) * NURBS_der000(:, j) -     sum003(j) * NURBS_der010(:, j) - 3 * sum012(j) * NURBS_der001(:, j) - 3 * sum002(j) * NURBS_der011(:, j) - 3 * sum011(j) * NURBS_der002(:, j) - 3 * sum001(j) * NURBS_der012(:, j) -     sum010(j) * NURBS_der003(:, j));
            NURBS_der004(:, j) = sum000_inv(j) * (Bspline_der004(:, j) - sum004(j) * NURBS_der000(:, j) - 4 * sum003(j) * NURBS_der001(:, j) - 6 * sum002(j) * NURBS_der002(:, j) - 4 * sum001(j) * NURBS_der003(:, j));
        end
        
        
    end
end


%{
function code = code_to_generate_coefficients(k1, k2, k3)
    k = k1 + k2 + k3;
    
    code = sprintf('NURBS_der%d%d%d(:, j) = sum000_inv(j) * (Bspline_der%d%d%d(:, j)', k1, k2, k3, k1, k2, k3);
    
    for j3 = 0 : min(k - 1, k3)
        for j2 = 0 : min(k - 1 - j3, k2)
            for j1 = 0 : min(k - 1 - j3 - j2, k1)
                coefficient = nchoosek(k1, j1) * nchoosek(k2, j2) * nchoosek(k3, j3);
                
                if (coefficient == 1)
                    if (j1 == 0 && j2 == 0 && j3 == 0)
                        code = sprintf('%s - sum%d%d%d(j) * NURBS_der%d%d%d(:, j)', code, k1 - j1, k2 - j2, k3 - j3, j1, j2, j3);
                    else
                        code = sprintf('%s -     sum%d%d%d(j) * NURBS_der%d%d%d(:, j)', code, k1 - j1, k2 - j2, k3 - j3, j1, j2, j3);
                    end
                else
                    code = sprintf('%s - %d * sum%d%d%d(j) * NURBS_der%d%d%d(:, j)', code, coefficient, k1 - j1, k2 - j2, k3 - j3, j1, j2, j3);
                end
            end
        end
    end
    
    code = sprintf('%s);', code);
end
%}