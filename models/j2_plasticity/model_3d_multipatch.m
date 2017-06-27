%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 3D J2 plasticity for isotropic and
%  homogeneous materials. Linear isotropic and kinematic hardening can be
%  accomodated.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_3d.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the directory to the global and patch
%          assembly files
%      path_to_results_directory is the directory to the results file
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement fields (.mat files)
%--------------------------------------------------------------------------
function model_3d_multipatch(path_to_assembly_directory, path_to_results_directory)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  3D J2 plasticity for isotropic, homogeneous materials.\n\n');
    
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches'              , ...
         'numNodesBeforePatch'     , ...
         'numMatrixEntries'        , ...
         'numDOFs'                 , ...
         'numDOFsPerNode'          , ...
         'GN_array'                , ...
         ...
         'material_sigmaY'         , ...
         'material_H'              , ...
         'isotropicHardeningFactor', ...
         'numLoadSteps'            , ...
         'numStressPoints'         , ...
         'numTimeStepsBetweenSaves', ...
         'maxNewtonsMethod'        , ...
         'tolNewtonsMethod');
    
    
    %----------------------------------------------------------------------
    %  Initialize the fields
    %  
    %  Note, we use Voigt notation and view the 3 x 3 strain and stress
    %  tensors as 6 x 1 vectors. The entries for the (total) strain and
    %  the plastic strain vectors are as follows:
    %  
    %    strain = [    strain_11;
    %                  strain_22;
    %                  strain_33;
    %              2 * strain_23;
    %              2 * strain_13;
    %              2 * strain_12];    (similarly for plastic strain)
    %  
    %  The entries of the stress vector are as follows:
    %  
    %    stress = [stress_11;
    %              stress_22;
    %              stress_33;
    %              stress_23;
    %              stress_13;
    %              stress_12];
    %  
    %----------------------------------------------------------------------
    % Assume zero displacements at t = 0
    u = zeros(numDOFs, 1);
    
    % Stress tensor at the next time step
    temp = zeros(6, numStressPoints);
    stress = temp;
    
    % Backstress tensor at the next time step
%    backstress_old = temp;
%    backstress     = temp;
    
    % Strain tensor at the next time step
    strain = temp;
    
    % Plastic strain tensors at the current and next time steps
    strain_pl_old = temp;
    strain_pl     = temp;
    clear temp;
    
    % Equivalent plastic strain at the current and next time steps
    strain_pl_eq_old = zeros(numStressPoints, 1);
    strain_pl_eq     = zeros(numStressPoints, 1);
    
    % Row, column, value arrays for the consistent tangent matrix
    temp = zeros(numMatrixEntries, 1);
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
            
    clear temp;
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over time (n)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for n = 1 : numTimeSteps
        fprintf('\n');
        fprintf('  Time index n = %d\n', n);
        
        % Load the BCs file
        load(sprintf('%sfile_bc_time%06.0f', path_to_assembly_directory, n), ...
             'ID_array', ...
             'BCU_array', ...
             'BCF_array', ...
             'numUnknownDOFs', ...
             'time_current');
        
        
        % Flag that indicates whether any of the material points is in
        % the plastic state
        isPlastic = 0;
        
        
        %------------------------------------------------------------------
        %  Set the BCs (provided by the user)
        %------------------------------------------------------------------
        % Displacement increment vector due to increments in displacements
        Delta_u = zeros(numDOFs, 1);
        Delta_u(BCU_array(:, 1)) = BCU_array(:, 2);
        
        % External force vector due to increments in body force and
        % tractions
        f_external = zeros(numDOFs, 1);
        f_external(BCF_array(:, 1)) = BCF_array(:, 2);
        
        % Indices for unknown displacements and unknown forces
        index_u = (1 : numUnknownDOFs)';
        index_f = ((numUnknownDOFs + 1) : numDOFs)';
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over Newton's method (k)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Initialize the new displacement vector
        u_new = u;
        
        % Flag that indicates whether we are close to convergence
        closeToConvergence = 0;
        
        % Norm of the residual vector
        residual_norm = 0;
        
        
        for k = 0 : maxNewtonsMethod
            % Initialize the internal force
            f_internal = zeros(numDOFs, 1);
            
            % Index that was last used to set the entry of K
            lastIndex_for_K = 0;
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over patches (p)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            for p = 1 : numPatches
                % Load the patch assembly file
                load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, p), ...
                     'elementType', ...
                     'material_E', ...
                     'material_nu', ...
                     'p1', ...
                     'p2', ...
                     'p3', ...
                     'nodes', ...
                     'numNodesPerElement', ...
                     'numDOFsPerElement', ...
                     'bezierExtractions1', ...
                     'bezierExtractions2', ...
                     'bezierExtractions3', ...
                     'numElements1', ...
                     'numElements2', ...
                     'numElements3', ...
                     'elementSizes1', ...
                     'elementSizes2', ...
                     'elementSizes3', ...
                     'IEN_array', ...
                     'numQuadraturePoints');
                
                switch (lower(elementType))
                    case {'bspline', 'b-spline'}
                        elementType = 0;
                        
                    case 'nurbs'
                        elementType = 1;
                end
                
                
                % Initialize the indices of the matrix entries
                index_for_K = ((lastIndex_for_K + 1) : (lastIndex_for_K + numDOFsPerElement))';
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Pre-process
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Some useful constants
                constant_p1p1 = p1 + 1;
                constant_p2p1 = p2 + 1;
                constant_p3p1 = p3 + 1;
                
                % Set the quadrature rule
                [z1, z2, z3, w] = set_3d_gauss_quadrature_for_bernstein(numQuadraturePoints);
                numQuadraturePointsPerElement = prod(numQuadraturePoints);
                
                % Evaluate the Bernstein polynomials for direction 1 at
                % the quadrature points
                Bernstein1_der0 = zeros(constant_p1p1, numQuadraturePoints(1));
                Bernstein1_der1 = zeros(constant_p1p1, numQuadraturePoints(1));
                
                for j = 1 : numQuadraturePoints(1)
                    temp = eval_1d_bernstein_der(z1(j), p1);
                    
                    Bernstein1_der0(:, j) = temp(:, 1);
                    Bernstein1_der1(:, j) = temp(:, 2);
                end
                
                % Evaluate the Bernstein polynomials for direction 2 at
                % the quadrature points
                Bernstein2_der0 = zeros(constant_p2p1, numQuadraturePoints(2));
                Bernstein2_der1 = zeros(constant_p2p1, numQuadraturePoints(2));
                
                for j = 1 : numQuadraturePoints(2)
                    temp = eval_1d_bernstein_der(z2(j), p2);
                    
                    Bernstein2_der0(:, j) = temp(:, 1);
                    Bernstein2_der1(:, j) = temp(:, 2);
                end
                
                % Evaluate the Bernstein polynomials for direction 3 at
                % the quadrature points
                Bernstein3_der0 = zeros(constant_p3p1, numQuadraturePoints(3));
                Bernstein3_der1 = zeros(constant_p3p1, numQuadraturePoints(3));
                
                for j = 1 : numQuadraturePoints(3)
                    temp = eval_1d_bernstein_der(z3(j), p3);
                    
                    Bernstein3_der0(:, j) = temp(:, 1);
                    Bernstein3_der1(:, j) = temp(:, 2);
                end
                
                clear z1 z2 z3 temp;
                
                
                %----------------------------------------------------------
                %  Build the elasticity matrix
                %----------------------------------------------------------
                % Bulk modulus
                material_K = material_E / (3 * (1 - 2 * material_nu));
                
                % Shear modulus
                material_G = material_E / (2 * (1 + material_nu));
                
                % D11, D22, D33 entries
                D11 = material_K + 4/3 * material_G;
                
                % D12, D13, D23, D21, D31, D32 entries
                D12 = material_K - 2/3 * material_G;
                
                % D44, D55, D66 entries
                D33 = material_G;
                
                % Set the elasticity matrix
                D_el = [D11, D12, D12,   0,   0,   0; ...
                        D12, D11, D12,   0,   0,   0; ...
                        D12, D12, D11,   0,   0,   0; ...
                          0,   0,   0, D33,   0,   0; ...
                          0,   0,   0,   0, D33,   0; ...
                          0,   0,   0,   0,   0, D33];
                
                
                %----------------------------------------------------------
                %  Build the LM array
                %----------------------------------------------------------
                LM_array = build_lm_array(IEN_array, ID_array, GN_array, numNodesBeforePatch(p));
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   End: Pre-process
                % ---------------------------------------------------------
                %----------------------------------------------------------
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over elements (e, e1, e2, e3)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Counters for the element and the quadrature point
                e = 1;
                q = 1;
                
                
                for e3 = 1 : numElements3
                    % Find the Bezier extraction matrix
                    bezierExtractions3_e = bezierExtractions3(:, :, e3);
                    
                    % Evaluate the map derivative dt/dzeta (constant)
                    dt_dzeta = 1 / elementSizes3(e3);
                    
                    % Evaluate the univariate B-splines
                    Bspline3_der0 = bezierExtractions3_e * Bernstein3_der0;
                    Bspline3_der1 = (dt_dzeta * bezierExtractions3_e) * Bernstein3_der1;
                    
                    
                    for e2 = 1 : numElements2
                        % Find the Bezier extraction matrix
                        bezierExtractions2_e = bezierExtractions2(:, :, e2);
                        
                        % Evaluate the map derivative dt/deta (constant)
                        dt_deta = 1 / elementSizes2(e2);
                        
                        % Evaluate the univariate B-splines
                        Bspline2_der0 = bezierExtractions2_e * Bernstein2_der0;
                        Bspline2_der1 = (dt_deta * bezierExtractions2_e) * Bernstein2_der1;
                        
                        
                        for e1 = 1 : numElements1
                            % Find the Bezier extraction matrix
                            bezierExtractions1_e = bezierExtractions1(:, :, e1);
                            
                            % Evaluate the map derivative dt/dxi (constant)
                            dt_dxi = 1 / elementSizes1(e1);
                            
                            % Evaluate the univariate B-splines
                            Bspline1_der0 = bezierExtractions1_e * Bernstein1_der0;
                            Bspline1_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
                            
                            
                            %----------------------------------------------
                            %  Evaluate the basis functions in the
                            %  parametric domain
                            %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
                            %----------------------------------------------
                            % For B-splines
                            if (elementType == 0)
                                % Find the positions of the nodes
                                nodes_e = nodes(IEN_array(:, e), :)';
                                
                                % Evaluate the B-splines
                                basis_der100 = kron(Bspline3_der0, kron(Bspline2_der0, Bspline1_der1));
                                basis_der010 = kron(Bspline3_der0, kron(Bspline2_der1, Bspline1_der0));
                                basis_der001 = kron(Bspline3_der1, kron(Bspline2_der0, Bspline1_der0));
                                
                            % For NURBS
                            elseif (elementType == 1)
                                % Find the positions and weights of the nodes
                                nodes_e = nodes(IEN_array(:, e), [1; 2; 3])';
                                w_e = nodes(IEN_array(:, e), 4);
                                
                                % Evaluate the B-splines
                                Bspline_der000 = kron(Bspline3_der0, kron(Bspline2_der0, Bspline1_der0));
                                Bspline_der100 = kron(Bspline3_der0, kron(Bspline2_der0, Bspline1_der1));
                                Bspline_der010 = kron(Bspline3_der0, kron(Bspline2_der1, Bspline1_der0));
                                Bspline_der001 = kron(Bspline3_der1, kron(Bspline2_der0, Bspline1_der0));
                                
                                % Evaluate the NURBS
                                [~, basis_der100, basis_der010, basis_der001] = eval_3d_nurbs_der1(w_e, Bspline_der000, Bspline_der100, Bspline_der010, Bspline_der001);
                                
                            else
                                fprintf('  Error: Element type must be B-splines or NURBS. The problem will be left unsolved.\n\n');
                                
                                quit;
                            end
                            
                            
                            %----------------------------------------------
                            %  Evaluate the map derivatives
                            %  (x = [x1; x2; x3])
                            %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                            %----------------------------------------------
                            dx_dxi = nodes_e * basis_der100;
                            dx_deta = nodes_e * basis_der010;
                            dx_dzeta = nodes_e * basis_der001;
                            
                            
                            
                            %----------------------------------------------
                            % ---------------------------------------------
                            %   Begin: Loop over quadrature points (q, q_e)
                            % ---------------------------------------------
                            %----------------------------------------------
                            % Initialize the element consistent tangent
                            % matrix and the element internal force vector
                            K_e = zeros(numDOFsPerElement);
                            f_e = zeros(numDOFsPerElement, 1);
                            
                            % Get the current Newton's guess for the
                            % solution vector
                            index_equation = LM_array(:, e);
                            u_e = u_new(index_equation);
                            
                            
                            for q_e = 1 : numQuadraturePointsPerElement
                                %------------------------------------------
                                %  Form the Jacobian matrix
                                %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                                %------------------------------------------
                                JacobianMatrix = [dx_dxi(:, q_e), ...
                                                  dx_deta(:, q_e), ...
                                                  dx_dzeta(:, q_e)]';
                                
                                % Evaluate the Jacobian
                                Jacobian = det(JacobianMatrix) / (dt_dxi * dt_deta * dt_dzeta);
                                
                                if (Jacobian <= 0)
                                    fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                                    
                                    quit;
                                end
                                
                                
                                %------------------------------------------
                                %  Evaluate the basis functions in the
                                %  physical domain
                                %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                                %------------------------------------------
                                basis_der1 = [basis_der100(:, q_e), ...
                                              basis_der010(:, q_e), ...
                                              basis_der001(:, q_e)]';
                                
                                basis_physical_der1 = JacobianMatrix \ basis_der1;
                                
                                
                                %------------------------------------------
                                %  Evaluate the element B matrix, which
                                %  gives the element strains when multiplied
                                %  by the element displacement vector
                                %------------------------------------------
                                % Initialize the B matrix
                                B = zeros(6, numDOFsPerElement);
                                
                                % Column indices for the B matrix
                                temp1 = 1;
                                temp2 = 2;
                                temp3 = 3;
                                
                                for a = 1 : numNodesPerElement
                                    B(1, temp1) = basis_physical_der1(1, a);
                                    B(5, temp1) = basis_physical_der1(3, a);
                                    B(6, temp1) = basis_physical_der1(2, a);
                                    
                                    B(2, temp2) = basis_physical_der1(2, a);
                                    B(4, temp2) = basis_physical_der1(3, a);
                                    B(6, temp2) = basis_physical_der1(1, a);
                                    
                                    B(3, temp3) = basis_physical_der1(3, a);
                                    B(4, temp3) = basis_physical_der1(2, a);
                                    B(5, temp3) = basis_physical_der1(1, a);
                                    
                                    temp1 = temp1 + numDOFsPerNode;
                                    temp2 = temp2 + numDOFsPerNode;
                                    temp3 = temp3 + numDOFsPerNode;
                                end
                                
                                
                                %------------------------------------------
                                %  Evaluate the trial strains and stresses
                                %------------------------------------------
                                % Evaluate the strain tensor
                                strain(:, q) = B * u_e;
                                
                                % Evaluate the elastic strain tensor and
                                % its trace
                                strain_el = strain(:, q) - strain_pl_old(:, q);
                                strain_el_tr = strain_el(1) + strain_el(2) + strain_el(3);
                                
                                % Evaluate the deviatoric elastic strain tensor
                                strain_el_dev = strain_el;
                                strain_el_dev(1) = strain_el_dev(1) - strain_el_tr / 3;
                                strain_el_dev(2) = strain_el_dev(2) - strain_el_tr / 3;
                                strain_el_dev(3) = strain_el_dev(3) - strain_el_tr / 3;
                                
                                % Evaluate the deviatoric stress tensor,
                                % relative to the backstress tensor
                                stress_dev = (2 * material_G) * strain_el_dev;
                                stress_dev(4) = stress_dev(4) / 2;
                                stress_dev(5) = stress_dev(5) / 2;
                                stress_dev(6) = stress_dev(6) / 2;
                                
                                
                                %------------------------------------------
                                %  Check the J2 yield condition
                                %------------------------------------------
                                % Evaluate the norm of the deviatoric
                                % stress tensor
                                stress_dev_norm = sqrt(stress_dev(1)^2 + stress_dev(2)^2 + stress_dev(3)^2 + 2 * (stress_dev(4)^2 + stress_dev(5)^2 + stress_dev(6)^2));
                                
                                % Evaluate the yield function
                                yield_function = stress_dev_norm - sqrt(2/3) * (material_sigmaY + material_H * strain_pl_eq_old(q));
                                
                                
                                %------------------------------------------
                                %  Elastic state
                                %------------------------------------------
                                if (yield_function <= 0)
                                    %--------------------------------------
                                    %  Accept the trial state
                                    %--------------------------------------
                                    % Do not update the plastic strain
                                    % tensor
                                    strain_pl(:, q) = strain_pl_old(:, q);
                                    
                                    % Do not update the equivalent plastic
                                    % strain
                                    strain_pl_eq(q) = strain_pl_eq_old(q);
                                    
                                    % Evaluate the stresses using the
                                    % elastic matrix
                                    stress(:, q) = D_el * strain_el;
                                    
                                    
                                    %--------------------------------------
                                    %  Evaluate the consistent tangent
                                    %  matrix and the internal force vector
                                    %--------------------------------------
                                    K_e = K_e + B' * (w(q_e) * Jacobian * D_el) * B;
                                    f_e = f_e + B' * (w(q_e) * Jacobian * stress(:, q));
                                    
                                    
                                %------------------------------------------
                                %  Plastic state
                                %------------------------------------------
                                else
                                    % Set the flag to 1
                                    isPlastic = 1;
                                    
                                    %--------------------------------------
                                    %  Radially return to the yield surface
                                    %--------------------------------------
                                    % Evaluate the consistency parameter
                                    Delta_gamma = yield_function / (2 * (material_G + material_H / 3));
                                    
                                    % Evaluate the normal, which is the unit
                                    % deviatoric stress tensor
                                    normal = stress_dev / stress_dev_norm;
                                    
                                    % Update the plastic strain
                                    strain_pl(:, q) = strain_pl_old(:, q) + Delta_gamma * [    normal(1); ...
                                                                                               normal(2); ...
                                                                                               normal(3); ...
                                                                                           2 * normal(4); ...
                                                                                           2 * normal(5); ...
                                                                                           2 * normal(6)];
                                    
                                    % Update the equivalent plastic strain
                                    strain_pl_eq(q) = strain_pl_eq_old(q) + sqrt(2/3) * Delta_gamma;
                                    
                                    % Update the stresses
                                    temp1 = 2 * material_G * Delta_gamma;
                                    
                                    stress(:, q) = stress_dev - temp1 * normal;
                                    stress(1, q) = stress(1, q) + material_K * strain_el_tr;
                                    stress(2, q) = stress(2, q) + material_K * strain_el_tr;
                                    stress(3, q) = stress(3, q) + material_K * strain_el_tr;
                                    
                                    
                                    %--------------------------------------
                                    %  Evaluate the consistent tangent
                                    %  matrix and the internal force vector
                                    %--------------------------------------
                                    % Here, the constants theta and theta_bar
                                    % are multiplied by (2 * material_G)
                                    temp2 = (2 * material_G) * (1 - temp1 / stress_dev_norm);
                                    temp3 = (2 * material_G) * (1 - temp2 / (2 * material_G) - material_G / (material_G + material_H / 3));
                                    
                                    % D11, D22, D33 entries
                                    D11 = material_K + 2 * temp2 / 3;
                                    
                                    % D12, D13, D23, D21, D31, D32 entries
                                    D12 = material_K - temp2 / 3;
                                    
                                    % D44, D55, D66 entries
                                    D33 = temp2 / 2;
                                    
                                    % Set the elastoplastic matrix
                                    D_ep = [D11, D12, D12,   0,   0,   0; ...
                                            D12, D11, D12,   0,   0,   0; ...
                                            D12, D12, D11,   0,   0,   0; ...
                                              0,   0,   0, D33,   0,   0; ...
                                              0,   0,   0,   0, D33,   0; ...
                                              0,   0,   0,   0,   0, D33] ...
                                         + temp3 * (normal * normal');
                                    
                                    K_e = K_e + B' * (w(q_e) * Jacobian * D_ep) * B;
                                    f_e = f_e + B' * (w(q_e) * Jacobian * stress(:, q));
                                    
                                    
                                end
                                
                                
                                % Increment the counter for quadrature point
                                q = q + 1;
                            end
                            
                            
                            %----------------------------------------------
                            % ---------------------------------------------
                            %   End: Loop over quadrature points
                            % ---------------------------------------------
                            %----------------------------------------------
                            
                            
                            %----------------------------------------------
                            %  Global assembly
                            %----------------------------------------------
                            % Add the element matrix
                            for j = 1 : numDOFsPerElement
                                rows_for_K(index_for_K)    = index_equation;
                                columns_for_K(index_for_K) = index_equation(j) * ones(numDOFsPerElement, 1);
                                values_for_K(index_for_K)  = K_e(:, j);
                                
                                index_for_K = index_for_K + numDOFsPerElement;
                            end
                            
                            % Add the element RHS vector
                            f_internal(index_equation) = f_internal(index_equation) + f_e;
                            
                            % Increment the counter for element
                            e = e + 1;
                        end
                    end
                end
                
                
                % Update the index last used to set the entry of K
                lastIndex_for_K = index_for_K(numDOFsPerElement);
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   End: Loop over elements
                % ---------------------------------------------------------
                %----------------------------------------------------------
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over patches
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            
            %--------------------------------------------------------------
            %  Check convergence
            %--------------------------------------------------------------
            if (k > 0)
                fprintf('\n');
                fprintf('  Close to convergence = %d\n', closeToConvergence);
                fprintf('  Relative error       = %.4e\n', relativeError);
                fprintf('  Residual norm        = %.4e\n', residual_norm);
                
                if (closeToConvergence == 1)
                    if (residual_norm <= residual_norm0)
                        fprintf('\n');
                        fprintf('  Plastic state = %d\n', isPlastic);
                        fprintf('  Newton''s method converged at iteration %d.\n\n', k);
                        
                        break;
                        
                    elseif (k == maxNewtonsMethod)
                        fprintf('\n');
                        fprintf('  Plastic state = %d\n', isPlastic);
                        fprintf('  Newton''s method did not converge after %d iterations.\n\n', k);
                        
                        quit;
                        
                    end
                    
                elseif (k == maxNewtonsMethod)
                    fprintf('\n');
                    fprintf('  Plastic state = %d\n', isPlastic);
                    fprintf('  Newton''s method did not converge after %d iterations.\n\n', k);
                    
                    quit;
                    
                end
            end
            
            
            %--------------------------------------------------------------
            %  Solve for the Newton's increment
            %--------------------------------------------------------------
            % Assemble the consistent tangent matrix
            K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs, numDOFs);
            
            % Save the current residual norm
            residual_norm0 = residual_norm;
            
            % Form the new residual vector (note the negative sign)
            residual = f_external - f_internal;
            residual_norm = norm(residual(index_u), 2);
            
            
            % Apply a row preconditioner
            precond = spdiags(1./max(K, [], 2), 0, numDOFs, numDOFs);
            K = precond * K;
            residual = precond * residual;
            clear precond;
            
            
            % Save the current Newton's increment
            Delta_u0 = Delta_u;
            
            % Solve for the new Newton's increment
            Delta_u(index_u) = K(index_u, index_u) \ (residual(index_u) - K(index_u, index_f) * Delta_u(index_f));
            u_new = u + Delta_u;
            
            
            % Compute the relative error in the Newton's increment
            denominator = norm(Delta_u, 2) + norm(Delta_u0, 2);
            relativeError = norm(Delta_u - Delta_u0, 2) / denominator;
            
            if (isPlastic == 0)
                tolerance = tolNewtonsMethod;
            else
                tolerance = 1e0;
            end
            
            if (relativeError < tolerance || denominator < 1e-15)
                closeToConvergence = 1;
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over Newton's method
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        
        %------------------------------------------------------------------
        %  Save the results
        %------------------------------------------------------------------
        u = u_new;
        strain_pl_old = strain_pl;
        strain_pl_eq_old = strain_pl_eq;
        
        if (mod(n, numTimeStepsBetweenSaves) == 0)
            save(sprintf('%sfile_results_time%06.0f', path_to_results_directory, n), 'time_current', 'u', 'stress', 'strain', 'strain_pl', 'strain_pl_eq', 'isPlastic', '-v7.3');
        end
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over time
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end
