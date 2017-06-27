%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 2D linear elasticity (plane strain)
%  for isotropic and homogeneous materials.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_2d_plane_strain.m output
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
%  1. Coefficients for the displacement fields (.mat file)
%--------------------------------------------------------------------------
function model_2d_plane_strain(path_to_assembly_directory, path_to_results_directory)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  2D linear elasticity (plane stress) for isotropic and homogeneous materials.\n\n');
    
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches', ...
         'numNodesBeforePatch', ...
         'numMatrixEntries', ...
         'numDOFs', ...
         'numDOFsPerNode', ...
         'GN_array', ...
         ...
         'material_thickness');
    
    % Load the BCs file
    load(sprintf('%sfile_bc', path_to_assembly_directory), ...
         'ID_array', ...
         'BCU_array', ...
         'BCF_array', ...
         'numUnknownDOFs');
    
    
    % Initialize the row, column, value arrays for the elasticity matrix
    temp = zeros(numMatrixEntries, 1);
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
    lastIndex_for_K = 0;
    
    clear temp;
    
    % Initialize the solution and RHS vectors
    u = zeros(numDOFs, 1);
    f = zeros(numDOFs, 1);
    
    % Prescribe the known values
    u(BCU_array(:, 1)) = BCU_array(:, 2);
    f(BCF_array(:, 1)) = BCF_array(:, 2);
    
    % Indices of unknown displacements and unknown forces
    index_u = (1 : numUnknownDOFs)';
    index_f = ((numUnknownDOFs + 1) : numDOFs)';
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over patches (p)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        % Load the patch assembly file
        load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, p), ...
             'description', ...
             'elementType', ...
             'material_E', ...
             'material_nu', ...
             'p1', ...
             'p2', ...
             'nodes', ...
             'numNodesPerElement', ...
             'numDOFsPerElement', ...
             'bezierExtractions1', ...
             'bezierExtractions2', ...
             'numElements1', ...
             'numElements2', ...
             'elementSizes1', ...
             'elementSizes2', ...
             'IEN_array', ...
             'numQuadraturePoints');
        
        fprintf('\n');
        fprintf('  Patch %d: %s\n', p, description);
        fprintf('    Degree of the %s basis functions: (%d, %d)\n', elementType, p1, p2);
        fprintf('    Number of elements in the parametric domain: (%d, %d)\n', numElements1, numElements2);
        
        switch (lower(elementType))
            case {'bspline', 'b-spline'}
                elementType = 0;
                
            case 'nurbs'
                elementType = 1;
                
        end
        
        
        % Initialize the indices of the matrix entries
        index_for_K = ((lastIndex_for_K + 1) : (lastIndex_for_K + numDOFsPerElement))';
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Pre-process
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Some useful constants
        constant_p1p1 = p1 + 1;
        constant_p2p1 = p2 + 1;
        
        % Set the quadrature rule
        [z1, z2, w] = set_2d_gauss_quadrature_for_bernstein(numQuadraturePoints);
        numQuadraturePointsPerElement = prod(numQuadraturePoints);
        
        % Evaluate the Bernstein polynomials for direction 1 at the
        % quadrature points
        Bernstein1_der0 = zeros(constant_p1p1, numQuadraturePoints(1));
        Bernstein1_der1 = zeros(constant_p1p1, numQuadraturePoints(1));
        
        for j = 1 : numQuadraturePoints(1)
            temp = eval_1d_bernstein_der(z1(j), p1);
            
            Bernstein1_der0(:, j) = temp(:, 1);
            Bernstein1_der1(:, j) = temp(:, 2);
        end
        
        % Evaluate the Bernstein polynomials for direction 2 at the
        % quadrature points
        Bernstein2_der0 = zeros(constant_p2p1, numQuadraturePoints(2));
        Bernstein2_der1 = zeros(constant_p2p1, numQuadraturePoints(2));
        
        for j = 1 : numQuadraturePoints(2)
            temp = eval_1d_bernstein_der(z2(j), p2);
            
            Bernstein2_der0(:, j) = temp(:, 1);
            Bernstein2_der1(:, j) = temp(:, 2);
        end
        
        clear z1 z2 temp;
        
        
        %------------------------------------------------------------------
        %  Build the elasticity matrix
        %------------------------------------------------------------------
        % Bulk modulus
        material_K = material_E / (3 * (1 - 2 * material_nu));
        
        % Shear modulus
        material_G = material_E / (2 * (1 + material_nu));
        
        % D11, D22 entries
        D11 = material_K + 4/3 * material_G;
        
        % D12, D21 entries
        D12 = material_K - 2/3 * material_G;
        
        % D33 entry
        D33 = material_G;
        
        % Set the elasticity matrix
        D_el = material_thickness * [D11, D12,   0; ...
                                     D12, D11,   0; ...
                                       0,   0, D33];
        
        
        %------------------------------------------------------------------
        %  Build the LM array
        %------------------------------------------------------------------
        LM_array = build_lm_array(IEN_array, ID_array, GN_array, numNodesBeforePatch(p));
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Pre-process
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over elements (e, e1, e2, e3)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Counter for the element
        e = 1;
        
        
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
                
                
                %----------------------------------------------------------
                %  Evaluate the basis functions in the parametric domain
                %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                % For B-splines
                if (elementType == 0)
                    % Find the positions of the nodes
                    nodes_e = nodes(IEN_array(:, e), :)';
                    
                    % Evaluate the B-splines
                    basis_der10 = kron(Bspline2_der0, Bspline1_der1);
                    basis_der01 = kron(Bspline2_der1, Bspline1_der0);
                    
                % For NURBS
                elseif (elementType == 1)
                    % Find the positions and weights of the nodes
                    nodes_e = nodes(IEN_array(:, e), [1; 2])';
                    w_e = nodes(IEN_array(:, e), 3);
                    
                    % Evaluate the B-splines
                    Bspline_der00 = kron(Bspline2_der0, Bspline1_der0);
                    Bspline_der10 = kron(Bspline2_der0, Bspline1_der1);
                    Bspline_der01 = kron(Bspline2_der1, Bspline1_der0);
                    
                    % Evaluate the NURBS
                    [~, basis_der10, basis_der01] = eval_2d_nurbs_der1(w_e, Bspline_der00, Bspline_der10, Bspline_der01);
                    
                else
                    fprintf('  Error: Element type must be B-splines or NURBS. The problem will be left unsolved.\n\n');
                    
                    quit;
                end
                
                
                %----------------------------------------------------------
                %  Evaluate the map derivatives (x = [x1; x2])
                %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                %----------------------------------------------------------
                dx_dxi = nodes_e * basis_der10;
                dx_deta = nodes_e * basis_der01;
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over quadrature points (q_e)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Initialize the element matrix and element RHS vector
                K_e = zeros(numDOFsPerElement);
%               f_e = zeros(numDOFsPerElement, 1);
                
                % Get the global equation indices
                index_equation = LM_array(:, e);
                
                
                for q_e = 1 : numQuadraturePointsPerElement
                    %------------------------------------------------------
                    %  Form the Jacobian matrix
                    %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                    %------------------------------------------------------
                    JacobianMatrix = [dx_dxi(:, q_e), ...
                                      dx_deta(:, q_e)]';
                    
                    % Evaluate the Jacobian
                    Jacobian = det(JacobianMatrix) / (dt_dxi * dt_deta);
                    
                    if (Jacobian <= 0)
                       fprintf('Error: Jacobian is not positive for e = %d, q_e = %d\n\n', e, q_e);
                       
                       quit;
                    end
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the basis functions in the physical domain
                    %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                    %------------------------------------------------------
                    basis_der1 = [basis_der10(:, q_e), ...
                                  basis_der01(:, q_e)]';
                    
                    basis_physical_der1 = JacobianMatrix \ basis_der1;
                    
                    
                    %------------------------------------------------------
                    %  Evaluate the element B matrix, which gives the
                    %  element strains when multiplied by the element
                    %  displacement vector
                    %------------------------------------------------------
                    % Initialize the B matrix
                    B = zeros(3, numDOFsPerElement);
                    
                    % Column indices for the B matrix
                    temp1 = 1;
                    temp2 = 2;
                    
                    for a = 1 : numNodesPerElement
                        B(1, temp1) = basis_physical_der1(1, a);
                        B(3, temp1) = basis_physical_der1(2, a);
                        
                        B(2, temp2) = basis_physical_der1(2, a);
                        B(3, temp2) = basis_physical_der1(1, a);
                        
                        temp1 = temp1 + numDOFsPerNode;
                        temp2 = temp2 + numDOFsPerNode;
                    end
                    
                    K_e = K_e + B' * ((w(q_e) * Jacobian) * D_el) * B;
                end
                
                
                %----------------------------------------------------------
                %  Global assembly
                %----------------------------------------------------------
                % Add the element matrix
                for j = 1 : numDOFsPerElement
                    rows_for_K(index_for_K)    = index_equation;
                    columns_for_K(index_for_K) = index_equation(j) * ones(numDOFsPerElement, 1);
                    values_for_K(index_for_K)  = K_e(:, j);
                    
                    index_for_K = index_for_K + numDOFsPerElement;
                end
                
                % Add the element RHS vector
%               f(index_equation) = f(index_equation) + f_e;
                
                % Increment the counter for element
                e = e + 1;
            end
        end
        
        
        % Update the index last used to set the entry of K
        lastIndex_for_K = index_for_K(numDOFsPerElement);
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over elements
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over patches
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    fprintf('\n');
    fprintf('  Solving for the displacements...\n\n');
    
    % Assemble the stiffness matrix
    K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs, numDOFs);
    clear rows_for_K columns_for_K values_for_K index_for_K K_e f_e;
    
    % Apply a row preconditioner
    precond = spdiags(1./max(K, [], 2), 0, numDOFs, numDOFs);
    K = precond * K;
    f = precond * f;
    clear precond;
    
    % Solve for the unknown coefficients
    u(index_u) = K(index_u, index_u) \ (f(index_u) - K(index_u, index_f) * u(index_f));
    clear K f;
    
    
    %----------------------------------------------------------------------
    %  Save the results
    %----------------------------------------------------------------------
    save(sprintf('%sfile_results', path_to_results_directory), 'u', '-v7.3');
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end
