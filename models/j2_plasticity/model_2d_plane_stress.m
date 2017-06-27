%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine solves the problem of 2D plane stress, J2 plasticity for
%  isotropic and homogeneous materials. Linear isotropic and kinematic
%  hardening can be accomodated.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_2d_plane_stress.m output
%  
%  Note that,
%  
%      path_to_assembly_directory is the path to the assembly files directory
%      path_to_results_directory is the path to the results directory
%      restartTime is the time step at which we start the simulation
%  
%  
%  Output:
%  
%  1. Coefficients for the displacement fields (.mat files)
%--------------------------------------------------------------------------
function model_2d_plane_stress(path_to_assembly_directory, path_to_results_directory, restartTime)
    % Feedback for user
    fprintf('\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
    fprintf('  2D plane stress, J2 plasticity for isotropic, homogeneous materials.\n\n');
    
    
    % Load the global assembly file
    load(sprintf('%sfile_assembly_global', path_to_assembly_directory), ...
         'numPatches'              , ...
         'numNodesBeforePatch'     , ...
         'numMatrixEntries'        , ...
         'numDOFs'                 , ...
         'numDOFsPerNode'          , ...
         'GN_array'                , ...
         ...
         'material_thickness'      , ...
         ...
         'material_sigmaY'         , ...
         'material_H'              , ...
         'isotropicHardeningFactor', ...
         'numLoadSteps'            , ...
         'numStressPoints'         , ...
         'numTimeStepsBetweenSaves', ...
         'maxNewtonsMethod'        , ...
         'tolNewtonsMethod');
    
    % Maximum number of iterations for finding the consistency parameter
    maxNewtonsMethod_Dg = 10;
    
    % Load the patch assembly file
    load(sprintf('%sfile_assembly_patch%d', path_to_assembly_directory, 1), ...
         'elementType'       , ...
         'material_E'        , ...
         'material_nu'       , ...
         'p1'                , ...
         'p2'                , ...
         'nodes'             , ...
         'numNodesPerElement', ...
         'numDOFsPerElement' , ...
         'bezierExtractions1', ...
         'bezierExtractions2', ...
         'numElements1'      , ...
         'numElements2'      , ...
         'elementSizes1'     , ...
         'elementSizes2'     , ...
         'IEN_array'         , ...
         'numQuadraturePoints');
    
    % Read the element type
    switch (lower(elementType))
        case {'bspline', 'b-spline'}
            elementType = 0;
            
        case 'nurbs'
            elementType = 1;
            
        otherwise
            fprintf('\n');
            fprintf('  Error: Element type must be B-splines or NURBS. The problem will be left unsolved.\n\n');
            
            quit;
    end
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Pre-process
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    if (restartTime == 0)
        %------------------------------------------------------------------
        %  Initialize loading
        %------------------------------------------------------------------
        % Set the load step (m) and time step (n) to the beginning
        m = 1;
        n = 1;
        
        % Initialize the total displacements and rotations applied
        displacementSoFar = [0; 0];
        rotationSoFar = 0;
        
        % Initialize the increment level
        incrementLevel = 0;
        
        % Initialize the number of successive convergences
        numSuccessiveConvergences = 0;
        
        % Flag that indicates whether the load step is finished
        flag_isLoaded = 0;
        
        
        %------------------------------------------------------------------
        %  Initialize the fields
        %  
        %  We use Voigt notation and view the 3 x 3 strain and stress
        %  tensors as 6 x 1 vectors. The entries for the (total) strain
        %  and the plastic strain vectors are as follows:
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
        %  Furthermore, we view the strain and stress tensors that satisfy
        %  the plane stress condition as 3 x 1 vectors.
        %  
        %    strain_2d = [    strain_11;
        %                     strain_22;
        %                 2 * strain_12]; (similarly for plastic strain)
        %    
        %    stress_2d = [stress_11;
        %                 stress_22;
        %                 stress_12];
        %  
        %------------------------------------------------------------------
        temp1 = zeros(6, numStressPoints);
        temp2 = zeros(3, numStressPoints);
        
        % Assume zero displacements at t = 0
        u = zeros(numDOFs, 1);
        
        % Stress tensor at the current time step
        stress    = temp1;
        stress_2d = temp2;
        
        % Backstress tensor at the current and previous time steps
        backstress_2d     = temp1;
        backstress_2d_old = temp2;
        
        % Strain tensor at the current time step
        strain    = temp1;
        strain_2d = temp2;
        
        % Plastic strain tensor at the current and previous time steps
        strain_pl        = temp1;
        strain_pl_2d     = temp2;
        strain_pl_2d_old = temp2;
        
        % Equivalent plastic strain at the current and previous time steps
        strain_pl_eq     = zeros(numStressPoints, 1);
        strain_pl_eq_old = zeros(numStressPoints, 1);
        
        clear temp1 temp2;
        
    else
        %------------------------------------------------------------------
        %  Initialize loading
        %------------------------------------------------------------------
        % Load the results from a previous time step
        load(sprintf('%sfile_results_time%06.0f', path_to_results_directory, restartTime), ...
             'u'                        , ...
             'stress'                   , ...
             'backstress_2d'            , ...
             'strain'                   , ...
             'strain_pl'                , ...
             'strain_pl_eq'             , ...
             ...
             'loadStep'                 , ...
             'displacementSoFar'        , ...
             'rotationSoFar'            , ...
             'incrementLevel'           , ...
             'numSuccessiveConvergences', ...
             'flag_isLoaded'            , ...
             '-v7.3');
        
        % Set the load step (m) and time step (n)
        m = loadStep;
        n = restartTime + 1;
        
        
        %------------------------------------------------------------------
        %  Initialize the fields
        %------------------------------------------------------------------
        % Read the backstress tensor, plastic strain tensor, and equivalent
        % plastic strain at the previous time step
        backstress_old   = backstress;
        strain_pl_old    = strain_pl;
        strain_pl_eq_old = strain_pl_eq;
        
    end
    
    
    %----------------------------------------------------------------------
    %  Initialize the row, column, and value arrays for the consistent
    %  tangent matrix
    %----------------------------------------------------------------------
    temp = zeros(numMatrixEntries, 1);
    
    rows_for_K    = temp;
    columns_for_K = temp;
    values_for_K  = temp;
            
    clear temp;
    
    
    %----------------------------------------------------------------------
    %  Set quadrature rule
    %----------------------------------------------------------------------
    % Some useful constants
    constant_p1p1 = p1 + 1;
    constant_p2p1 = p2 + 1;
    
    % Set the quadrature rule
    [z1, z2, w] = set_2d_gauss_quadrature_for_bernstein(numQuadraturePoints);
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
    
    clear z1 z2 temp;
    
    
    %----------------------------------------------------------------------
    %  Build the elasticity matrix
    %----------------------------------------------------------------------
    % Bulk modulus
    material_K = material_E / (3 * (1 - 2 * material_nu));
    
    % Shear modulus
    material_G = material_E / (2 * (1 + material_nu));
    
    % D11, D22 entries
    D11 = material_E / (1 - material_nu^2);
    
    % D12, D21 entries
    D12 = material_E * material_nu / (1 - material_nu^2);
    
    % D33 entry
    D33 = material_G;
    
    % Set the elasticity matrix
    D_el = material_thickness * [D11, D12,   0; ...
                                 D12, D11,   0; ...
                                   0,   0, D33];
    
    
    %----------------------------------------------------------------------
    %  Build the inverse of the elasticity matrix
    %----------------------------------------------------------------------
    % C11, C22 entries
    C11 = 1 / material_E;
    
    % C12, C21 entries
    C12 = -material_nu / material_E;
    
    % C33 entry
    C33 = 2 * (1 + material_nu) / material_E;
    
    D_el_inv = 1/material_thickness * [C11, C12,   0; ...
                                       C12, C11,   0; ...
                                         0,   0, C33];
    
    
    %----------------------------------------------------------------------
    %  Build the P matrix, which maps the stress vector (satisfying plane
    %  stress condition) to the deviatoric stress vector
    %----------------------------------------------------------------------
    P = [ 2/3, -1/3, 0; ...
         -1/3,  2/3, 0; ...
            0,    0, 2];
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Pre-process
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   Begin: Loop over load steps (m)
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    % Some useful constants for plasticity
    constant_sq2d3 = sqrt(2/3);
    constant_2G = 2 * material_G;
    
    
    % Initialization for solution and RHS vectors
    globalVector  = zeros(numDOFs, 1);
    
    % Initialization for element matrix K_e
    elementMatrix = zeros(numDOFsPerElement);
    
    % Initialization for element vector f_e
    elementVector = zeros(numDOFsPerElement, 1);
    
    % Number of matrix entries that we compute for each element
    numMatrixEntriesPerElement = numDOFsPerElement^2;
    
    
    % Vector of ones for global assembly
    vector_of_ones = ones(numDOFsPerElement, 1);
    
    % Vector of indices for element matrices
    indices_for_elementMatrix = (1 : numMatrixEntriesPerElement)';
    
    
    for loadStep = m : numLoadSteps
        % Find the positions of the nodes at the beginning of the load step
        % (TODO)
        nodesAtBeginning = nodes(:, [1; 2])';
        
        
        % Load the load step file
        load(sprintf('%sfile_loadstep%d', path_to_assembly_directory, loadStep), ...
             'numBoundaries', ...
             'numTimeSteps');
        
        % Initialize the increment
        increment = (1 / numTimeSteps) * 2^incrementLevel;
        
        % Initialize the ID array
        ID_array = [];
        
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   Begin: Loop over time steps (n)
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        % Continue to loop until we reach the total displacement or the
        % total rotation that is specified by the load step
        while (flag_isLoaded ~= 1)
            fprintf('\n');
            fprintf('- Load step = %d, time step = %d\n', loadStep, n);
            fprintf('  Increment level = %d\n', incrementLevel);
            
            
            %--------------------------------------------------------------
            %  Set the displacement and force BCs
            %--------------------------------------------------------------
            BCs_displacement = [];
            BCs_force = [];
            
            for boundary = 1 : numBoundaries
                % Load the BC file
                load(sprintf('%sfile_loadstep%d_boundary%d', path_to_assembly_directory, loadStep, boundary), ...
                     'loadType', ...
                     'nodesOnBoundary', ...
                     'dofIndex', ...
                     'displacementTotal', ...
                     'rotationTotal', ...
                     'rotationAxis');
                
                numNodesOnBoundary = size(nodesOnBoundary, 1);
                numDOFIndex = size(dofIndex, 1);
                
                % Displacement only
                if (loadType == 1)
                    % Evaluate the incremental displacements for the nodes
                    % on the boundary
                    bcValue = increment * displacementTotal;
                    
                    BCs_displacement = [BCs_displacement; ...
                                        [kron(nodesOnBoundary, ones(numDOFIndex, 1)), repmat(dofIndex, numNodesOnBoundary, 1), repmat(bcValue, numNodesOnBoundary, 1)]];
                    
                % Rotation only
                elseif (loadType == 2)
                    % Evaluate the rotation matrix that maps the nodes on 
                    % the boundary at the beginning of the load step to
                    % those at the current time step
                    theta0 = rotationSoFar;
                    
                    constant_costheta = cos(theta0);
                    constant_sintheta = sin(theta0);
                    
                    rotationMatrix0 = [constant_costheta, -constant_sintheta; ...
                                       constant_sintheta,  constant_costheta];
                    
                    % Evaluate the rotation matrix that maps the nodes on 
                    % the boundary at the beginning of the load step to
                    % those at the next time step
                    theta1 = rotationSoFar + increment * rotationTotal;
                    
                    constant_costheta = cos(theta1);
                    constant_sintheta = sin(theta1);
                    
                    rotationMatrix0 = [constant_costheta, -constant_sintheta; ...
                                       constant_sintheta,  constant_costheta];
                    
                    
                    % Evaluate the incremental displacements for the nodes
                    % on the boundary
                    bcValue = (rotationMatrix1 - rotationMatrix0) * nodesAtBeginning(dofIndex, nodesOnBoundary);
                    
                    BCs_displacement = [BCs_displacement; ...
                                        [kron(nodesOnBoundary, ones(numDOFIndex, 1)), repmat(dofIndex, numNodesOnBoundary, 1), reshape(bcValue, numNodesOnBoundary * numDOFIndex, 1)]];
                    
                end
            end
            
            
            %--------------------------------------------------------------
            %  Build the ID, LM, and BC arrays
            %--------------------------------------------------------------
            % Unless the nodes on the boundaries change within a load step,
            % we only need to build the ID and LM arrays once
            if (isempty(ID_array))
                ID_array = build_id_array(BCs_displacement, numDOFsPerNode, GN_array);
                LM_array = build_lm_array(IEN_array, ID_array, GN_array, 0);
            end
            
            % We build the BC arrays at every time step since the BC value
            % differ from one time step to another. The equation indices
            % will stay the same, however.
            [BCU_array, BCF_array] = build_bc_array(BCs_displacement, BCs_force, ID_array);
            
            % Number of unknown degrees of freedom
            numUnknownDOFs = numDOFs - size(BCU_array, 1);
            
            % Indices for unknown displacements
            index_u = (1 : numUnknownDOFs)';
            index_f = ((numUnknownDOFs + 1) : numDOFs)';
            
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   Begin: Loop over Newton's method (k)
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            % By definition, the initial Newton's guess is the solution
            % vector from the previous time step
            u_new = u;
            
            % External force vector due to increments in body force and
            % tractions
            f_external = globalVector;
            f_external(BCF_array(:, 1)) = BCF_array(:, 2);
            
            
            % Initialize the Newton's increment
            u_increment = globalVector;
            
            
            % Flag that indicates whether any of the material points is in
            % the plastic state
            flag_isPlastic = 0;
            
            % Flag that indicates whether we are close to convergence
            flag_isCloseToConvergence = 0;
            
            % Norm of the residual vector
            residual_norm = 0;
            
            
            for k = 0 : maxNewtonsMethod
                % Initialize the internal force
                f_internal = globalVector;
                
                % Index that was last used to set the entry of K
                lastIndex_for_K = 0;
                
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   Begin: Loop over elements (e, e1, e2)
                % ---------------------------------------------------------
                %----------------------------------------------------------
                % Counters for the element and the quadrature point
                e = 1;
                q = 1;
                
                
                for e2 = 1 : numElements2
                    % Find the Bezier extraction matrix
                    bezierExtractions2_e = bezierExtractions2(:, :, e2);
                    
                    % Evaluate the map derivative dt/deta (constant)
                    dt_deta = 1 / elementSizes2(e2);
                    
                    % Evaluate the univariate B-splines
                    Bspline2_der0 =            bezierExtractions2_e  * Bernstein2_der0;
                    Bspline2_der1 = (dt_deta * bezierExtractions2_e) * Bernstein2_der1;
                    
                    
                    for e1 = 1 : numElements1
                        % Find the Bezier extraction matrix
                        bezierExtractions1_e = bezierExtractions1(:, :, e1);
                        
                        % Evaluate the map derivative dt/dxi (constant)
                        dt_dxi = 1 / elementSizes1(e1);
                        
                        % Evaluate the univariate B-splines
                        Bspline1_der0 =           bezierExtractions1_e  * Bernstein1_der0;
                        Bspline1_der1 = (dt_dxi * bezierExtractions1_e) * Bernstein1_der1;
                        
                        
                        %--------------------------------------------------
                        %  Evaluate the basis functions in the parametric domain
                        %  Matrix: (numNodesPerElement) x (numQuadraturePointsPerElement)
                        %--------------------------------------------------
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
                            w_e     = nodes(IEN_array(:, e), 3);
                            
                            % Evaluate the B-splines
                            Bspline_der00 = kron(Bspline2_der0, Bspline1_der0);
                            Bspline_der10 = kron(Bspline2_der0, Bspline1_der1);
                            Bspline_der01 = kron(Bspline2_der1, Bspline1_der0);
                            
                            % Evaluate the NURBS
                            [~, basis_der10, basis_der01] = eval_2d_nurbs_der1(w_e, Bspline_der00, Bspline_der10, Bspline_der01);
                            
                        end
                        
                        
                        %--------------------------------------------------
                        %  Evaluate the map derivatives (x = [x1; x2])
                        %  Matrix: (numDOFsPerNode) x (numQuadraturePointsPerElement)
                        %--------------------------------------------------
                        dx_dxi  = nodes_e * basis_der10;
                        dx_deta = nodes_e * basis_der01;
                        
                        
                        
                        %--------------------------------------------------
                        % -------------------------------------------------
                        %   Begin: Loop over quadrature points (q, q_e)
                        % -------------------------------------------------
                        %--------------------------------------------------
                        % Initialize the element consistent tangent matrix and
                        % the element internal force vector
                        K_e = elementMatrix;
                        f_e = elementVector;
                        
                        % Get the current Newton's guess for the incremental displacements
                        index_equation = LM_array(:, e);
                        u_e = u_new(index_equation);
                        
                        
                        for q_e = 1 : numQuadraturePointsPerElement
                            %----------------------------------------------
                            %  Form the Jacobian matrix
                            %  Matrix: (numDOFsPerNode) x (numDOFsPerNode)
                            %----------------------------------------------
                            JacobianMatrix = [dx_dxi(:, q_e), dx_deta(:, q_e)]';
                            
                            % Evaluate the Jacobian
                            Jacobian = det(JacobianMatrix) / (dt_dxi * dt_deta);
                            
                            if (Jacobian <= 0)
                                fprintf('\n');
                                fprintf('  Error: Jacobian is not positive for e = %d, q_e = %d. The problem will be left unsolved.\n\n', e, q_e);
                                
                                quit;
                            end
                            
                            
                            %----------------------------------------------
                            %  Evaluate the basis functions in the physical domain
                            %  Matrix: (numDOFsPerNode) x (numNodesPerElement)
                            %----------------------------------------------
                            basis_der1 = [basis_der10(:, q_e), basis_der01(:, q_e)]';
                            
                            basis_physical_der1 = JacobianMatrix \ basis_der1;
                            
                            
                            %----------------------------------------------
                            %  Evaluate the element B matrix, which gives the element
                            %  strains when multiplied by the element displacement vector
                            %----------------------------------------------
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
                                
                                temp1 = temp1 + 2;
                                temp2 = temp2 + 2;
                            end
                            
                            
                            %----------------------------------------------
                            %  Evaluate the trial strains and stresses
                            %----------------------------------------------
                            % Evaluate the strain tensor
                            strain_2d(:, q) = B * u_e;
                            
                            % Evaluate the elastic strain tensor
                            strain_el_2d = strain_2d(:, q) - strain_pl_2d_old(:, q);
                            
                            % Evaluate the stress tensor
                            stress_2d(:, q) = D_el * strain_el_2d;
                            
                            % Evaluate the stress tensor, relative to the backstress tensor
                            % (we implicitly map the relative stress to its deviator below)
                            stress_rel_2d = stress_2d(:, q) - backstress_2d_old(:, q);
                            
                            
                            %----------------------------------------------
                            %  Check the J2 yield condition
                            %----------------------------------------------
                            % Some useful constants for radial return
                            numerator1 = (stress_rel_2d(1) + stress_rel_2d(2))^2 / 6;
                            numerator2 = (stress_rel_2d(2) - stress_rel_2d(1))^2 / 2 + 2 * stress_rel_2d(3)^2;
                            
                            % Evaluate the norm of the deviatoric stress tensor
                            % Note that,
                            %   s : s = sigma : s
                            %         = stress_rel_2d' * P * stress_rel_2d
                            %         = (stress_rel_2d(1) + stress_rel_2d(2))^2 / 6 + (stress_rel_2d(2) - stress_rel_2d(1))^2 / 2 + 2 * stress_rel_2d(3)^2
                            stress_dev_norm = sqrt(numerator1 + numerator2);
                            
                            % Evaluate the hardening function
                            hardening_function = material_sigmaY + isotropicHardeningFactor * material_H * strain_pl_eq_old(q);
                            
                            % Evaluate the yield function
                            yield_function = stress_dev_norm^2 / 2 - hardening_function^2 / 3;
                            
                            
                            %----------------------------------------------
                            %  Elastic state
                            %----------------------------------------------
                            if (yield_function <= 0)
                                %------------------------------------------
                                %  Accept the trial state
                                %------------------------------------------
                                % Do not update the backstress tensor
                                if (isotropicHardeningFactor < 1)
                                    backstress_2d(:, q) = backstress_2d_old(:, q);
                                end
                                
                                % Do not update the plastic strain tensor
                                strain_pl_2d(:, q) = strain_pl_2d_old(:, q);
                                
                                % Do not update the equivalent plastic strain
                                strain_pl_eq(q) = strain_pl_eq_old(q);
                                
                                
                                %------------------------------------------
                                %  Evaluate the consistent tangent matrix and
                                %  the internal force vector
                                %------------------------------------------
                                K_e = K_e + B' * (w(q_e) * Jacobian * D_el) * B;
                                f_e = f_e + B' * (w(q_e) * Jacobian * stress_2d(:, q));
                                
                                
                            %----------------------------------------------
                            %  Plastic state
                            %----------------------------------------------
                            else
                                % Set the flag to 1
                                flag_isPlastic = 1;
                                
                                
                                %------------------------------------------
                                %  Evaluate the consistency parameter
                                %------------------------------------------
                                % Some useful constants for radial return
                                constant1 = 2/3 * (1 - isotropicHardeningFactor) * material_H;
                                constant2 = constant1 + material_E / (3 * (1 - material_nu));
                                constant3 = constant1 + constant_2G;
                                
                                
                                % Initial the Newton's guess for the consistency parameter (Dg)
                                Delta_gamma = 0;
                                
                                % Flag that indicates whether we are close to convergence
                                flag_closeToConvergence_Dg = 0;
                                
                                
                                for k_Dg = 0 : maxNewtonsMethod_Dg
                                    % Some useful constants for radial return
                                    denominator1 = 1 + Delta_gamma * constant2;
                                    denominator2 = 1 + Delta_gamma * constant3;
                                    
                                    
                                    % Evaluate the norm of the deviatoric stress tensor
                                    stress_dev_norm = sqrt(numerator1 / denominator1^2 + numerator2 / denominator2^2);
                                    
                                    % Evaluate the derivative of the square of the norm of the 
                                    % deviatoric stress tensor with respect to the consistency parameter
                                    stress_dev_norm_sq_der = (-2 * numerator1 / denominator1^3) * constant2 + (-2 * numerator2 / denominator2^3) * constant3;
                                    
                                    
                                    % Evaluate the hardening function
                                    hardening_function = material_sigmaY + isotropicHardeningFactor * material_H * (strain_pl_eq_old(q) + Delta_gamma * constant_sq2d3 * stress_dev_norm);
                                    
                                    % Evaluate the derivative of the square of the hardening function
                                    % with respect to the consistency parameter
                                    hardening_function_sq_der = (2 * constant_sq2d3 * hardening_function) * (isotropicHardeningFactor * material_H) * (stress_dev_norm + Delta_gamma * stress_dev_norm_sq_der / (2 * stress_dev_norm));
                                    
                                    
                                    %--------------------------------------
                                    %  Check for convergence
                                    %--------------------------------------
                                    if (k_Dg > 0)
                                        fprintf('\n');
                                        fprintf('  Residual = %.4e\n', residual_Dg_norm);
                                        
                                        if (flag_closeToConvergence_Dg == 1 && residual_Dg_norm < residual_Dg_norm0)
                                            fprintf('\n');
                                            fprintf('  Newton''s method for consistency parameter converged at iteration %d.\n\n', k_Dg);
                                            
                                            break;
                                            
                                        elseif (k_Dg == maxNewtonsMethod_Dg)
                                            fprintf('\n');
                                            fprintf('  Newton''s method for consistency parameter did not converge after %d iterations.\n\n', k_Dg);
                                            
                                            quit;
                                        end
                                    end
                                    
                                    
                                    %--------------------------------------
                                    %  Solve for the Newton's increment
                                    %--------------------------------------
                                    % Assemble the tangent matrix (the derivative of the residual)
                                    K_Dg = stress_dev_norm_sq_der / 2 - hardening_function_sq_der / 3;
                                    
                                    
                                    % Save the norm of the residual vector at the previous iteration
                                    residual_Dg_norm0 = residual_Dg_norm;
                                    
                                    % Form the residual vector
                                    residual_Dg = stress_dev_norm^2 / 2 - hardening_function^2 / 3;
                                    
                                    % Evaluate the norm of the residual vector at the current iteration
                                    residual_Dg_norm = abs(residual_Dg);
                                    
                                    
                                    % Solve for the Newton's increment
                                    Delta_gamma_increment = -residual_Dg / K_Dg;
                                    
                                    % Update the Newton's guess
                                    Delta_gamma = Delta_gamma + Delta_gamma_increment;
                                    
                                    
                                    % Evaluate the norm of the Newton's increment and
                                    % check if we are close to convergence
                                    increment_Dg_norm = abs(Delta_gamma_increment);
                                    
                                    if (increment_Dg_norm < 1e-6)
                                        flag_closeToConvergence_Dg = 1;
                                    end
                                end
                                
                                
                                %------------------------------------------
                                %  Build the modified elastic tangent matrix
                                %------------------------------------------
                                % Some useful constant for radial return
                                constant4 = 1 + Delta_gamma * constant1;
                                
                                D_modified = (D_el_inv + (Delta_gamma / constant4) * P) \ eye(3);
                                
                                
                                %------------------------------------------
                                %  Radially return to the yield surface
                                %------------------------------------------
                                % Update the stress tensor, relative to the backstress tensor
                                stress_rel_2d = (1 / constant4) * (D_modified * D_el_inv * stress_rel_2d);
                                
                                % Update the backstress tensor
                                if (isotropicHardeningFactor < 1)
                                    backstress_2d(:, q) = backstress_2d_old(:, q) + (constant4 - 1) * stress_rel_2d;
                                end
                                
                                % Update the stress tensor
                                stress_2d(:, q) = stress_rel_2d + backstress_2d(:, q);
                                
                                % Map the relative stress to its deviator
                                stress_rel_dev = P * stress_rel_2d;
                                
                                % Update the plastic strain tensor
                                strain_pl_2d(:, q) = strain_pl_2d_old(:, q) + Delta_gamma * stress_rel_dev;
                                
                                % Update the equivalent plastic strain
                                strain_pl_eq(q) = strain_pl_eq_old(q) + Delta_gamma * constant_sq2d3 * stress_dev_norm;
                                
                                
                                %------------------------------------------
                                %  Evaluate the consistent tangent matrix and
                                %  the internal force vector
                                %------------------------------------------
                                % Some useful constants for the consistent tangent matrix
                                constant5 = 1 - Delta_gamma * (2/3 * isotropicHardeningFactor * material_H);
                                constant6 = 2/3 * (constant4 / constant5) * material_H;
                                
                                % A misnomer, but we call this the normal tensor
                                normal = D_modified * stress_rel_dev;
                                
                                % Set the elastoplastic matrix
                                D_ep = D_modified - (normal * normal') / (stress_rel_dev' * (normal + constant6 * stress_rel_2d));
                                
                                K_e = K_e + B' * (w(q_e) * Jacobian * D_ep) * B;
                                f_e = f_e + B' * (w(q_e) * Jacobian * stress(:, q));
                                
                                
                            end
                            
                            
                            % Increment the counter for quadrature point
                            q = q + 1;
                        end
                        
                        
                        %--------------------------------------------------
                        % -------------------------------------------------
                        %   End: Loop over quadrature points
                        % -------------------------------------------------
                        %--------------------------------------------------
                        
                        
                        %--------------------------------------------------
                        %  Global assembly
                        %--------------------------------------------------
                        % Add the element matrix
                        index = lastIndex_for_K + indices_for_elementMatrix;
                        
                        rows_for_K(index)    = kron(vector_of_ones, index_equation);
                        columns_for_K(index) = kron(index_equation, vector_of_ones);
                        values_for_K(index)  = reshape(K_e, numMatrixEntriesPerElement, 1);
                        
                        index = index + numMatrixEntriesPerElement;
                        
                        
                        % Add the element RHS vector
                        f_internal(index_equation) = f_internal(index_equation) + f_e;
                        
                        
                        % Increment the counter for element
                        e = e + 1;
                    end
                end
                
                
                %----------------------------------------------------------
                % ---------------------------------------------------------
                %   End: Loop over elements
                % ---------------------------------------------------------
                %----------------------------------------------------------
                
                
                
                %----------------------------------------------------------
                %  Check for convergence
                %  
                %  We accept that the Newton's guess corresponds to an
                %  equilibrium state via a two-step check.
                %----------------------------------------------------------
                if (k > 0)
                    fprintf('\n');
                    fprintf('  Close to convergence = %d\n', flag_isCloseToConvergence);
                    fprintf('  Newton''s increment abs., 2-norm = %.4e\n', increment_norm);
                    fprintf('  Residual vector    abs., 2-norm = %.4e\n', residual_norm);
                    
                    % Check whether the current Newton's increment is small,
                    % and if so, whether the current residual is smaller than
                    % the previous residual
                    if (flag_isCloseToConvergence == 1 && residual_norm < residual_norm0)
                        fprintf('\n');
                        fprintf('  Newton''s method converged at iteration %d.\n', k);
                        fprintf('  Plastic state = %d\n\n', flag_isPlastic);
                        
                        % Set the flag to 1
                        flag_hasConverged = 1;
                        
                        % Record the number of successive convergences
                        if (incrementLevel < 0)
                            numSuccessiveConvergences = numSuccessiveConvergences + 1;
                        end
                        
                        break;
                        
                    elseif (k == maxNewtonsMethod)
                        fprintf('\n');
                        fprintf('  Newton''s method did not converge after %d iterations.\n', k);
                        fprintf('  Plastic state = %d\n\n', flag_isPlastic);
                        
                        % Set the flag to 0
                        flag_hasConverged = 0;
                        
                        % Reset the number of successive convergences
                        numSuccessiveConvergences = 0;
                        
                        break;
                        
                    end
                end
                
                
                %----------------------------------------------------------
                %  Solve for the Newton's increment
                %----------------------------------------------------------
                % Assemble the consistent tangent matrix
                K = sparse(rows_for_K, columns_for_K, values_for_K, numDOFs, numDOFs);                
                
                
                % Save the residual norm at the previous iteration
                residual_norm0 = residual_norm;
                
                % Form the residual vector (note the negative sign)
                residual = f_external - f_internal;
                
                % Evaluate the residual norm at the current iteration
                residual_norm = norm(residual(index_u), 2);
                
                
                % Apply a row preconditioner
                precond = spdiags(1./max(K, [], 2), 0, numDOFs, numDOFs);
                K = precond * K;
                residual = precond * residual;
                
                
                % At the first iteration, we tell the Newton's increment to
                % satisfy any incremental displacements on the boundary
                if (k == 0)
                    u_increment(BCU_array(:, 1)) = BCU_array(:, 2);
                    u_increment(index_u) = K(index_u, index_u) \ (residual(index_u) - K(index_u, index_f) * u_increment(index_f));
                    
                % At subsequent iterations, the Newton's guess already satisfies
                % the displacement BCs. Hence, we tell the Newton's increment to
                % stop imposing the incremental displacements on the bodundary
                else
                    if (k == 1)
                        u_increment(BCU_array(:, 1)) = 0;
                    end
                    u_increment(index_u) = K(index_u, index_u) \ residual(index_u);
                    
                end
                
                % Update the Newton's guess
                u_new = u_new + u_increment;
                
                
                % Evaluate the norm of the Newton's increment and check if
                % we are close to convergence
                increment_norm = norm(u_increment(index_u), 2);
                
                if (increment_norm < tolNewtonsMethod)
                    flag_isCloseToConvergence = 1;
                end
            end
            
            
            %--------------------------------------------------------------
            % -------------------------------------------------------------
            %   End: Loop over Newton's method
            % -------------------------------------------------------------
            %--------------------------------------------------------------
            
            
            %--------------------------------------------------------------
            %  Decide what to do next
            %--------------------------------------------------------------
            if (flag_hasConverged == 1)
                %----------------------------------------------------------
                %  Evaluate the out-of-plane stress and strain components
                %----------------------------------------------------------
                for q = 1 : numStressPoints
                    % Stress tensor
                    stress(:, q) = [stress_2d(1, q); ...
                                    stress_2d(2, q); ...
                                                  0; ...
                                                  0; ...
                                                  0; ...
                                    stress_2d(3, q)];
                    
                    % Strain tensor
                    strain(:, q) = [strain_2d(1, q); ...
                                    strain_2d(2, q); ...
                                    -material_nu / (1 - material_nu) * (strain_2d(1, q) + strain_2d(2, q)); ...
                                                  0; ...
                                                  0; ...
                                    strain_2d(3, q)];
                    
                    % Plastic strain tensor
                    strain_pl(:, q) = [strain_pl_2d(1, q); ...
                                       strain_pl_2d(2, q); ...
                                       -(strain_pl_2d(1, q) + strain_pl_2d(2, q)); ...
                                                        0; ...
                                                        0; ...
                                       strain_pl_2d(3, q)];
                    
                end
                
                
                if (loadType == 1)
                    % Update the total displacement applied
                    displacementSoFar = displacementSoFar + increment * displacementTotal;
                    
                    % Check if we are close to finishing the load step
                    if (displacementSoFar(1) >= displacementTotal(1) && displacementSoFar(2) >= displacementTotal(2))
                        flag_isLoaded = 1;
                    end
                    
                elseif (loadType == 2)
                    % Update the total rotation applied
                    rotationSoFar = rotationSoFar + increment * rotationTotal;
                    
                    % Check if we are close to finishing the load step
                    if (rotationSoFar >= rotationTotal)
                        flag_isLoaded = 1;
                    end
                    
                end
                
                
                % Save the results
                if (mod(n, numTimeStepsBetweenSaves) == 0 || flag_isLoaded == 1)
                    save(sprintf('%sfile_results_time%06.0f', path_to_results_directory, n), ...
                         'ID_array', ...
                         'u', 'stress', 'backstress_2d', 'strain', 'strain_pl', 'strain_pl_eq', ...
                         'loadStep', 'displacementSoFar', 'rotationSoFar', 'incrementLevel', 'numSuccessiveConvergences', 'flag_isLoaded', ...
                         '-v7.3');
                end
                
                
                % Update the time step
                n = n + 1;
                
                % Update the fields
                u = u_new;
                backstress_2d_old = backstress_2d;
                strain_pl_2d_old = strain_pl_2d;
                strain_pl_eq_old = strain_pl_eq;
                
                % Check whether to relax the increment
                if (incrementLevel < 0 && numSuccessiveConvergences == 4)
                    increment = 2 * increment;
                    incrementLevel = incrementLevel + 1;
                end
                
            else
                % Lower the increment
                increment = 0.5 * increment;
                incrementLevel = incrementLevel - 1;
                
                % If the increment level is deemed too low, then terminate
                % and consider the problem unsolved
                if (incrementLevel < -7)
                    fprintf('\n');
                    fprintf('  Error: Equilibrium could not be reached with the lowest increment. The problem will be left unsolved.\n\n');
                    
                    quit;
                end
                
            end
        end
        
        
        %------------------------------------------------------------------
        % -----------------------------------------------------------------
        %   End: Loop over time step
        % -----------------------------------------------------------------
        %------------------------------------------------------------------
        
        
        % Reset the values for the next load step
        displacementSoFar = [0; 0];
        rotationSoFar = 0;
        incrementLevel = 0;
        numSuccessiveConvergences = 0;
        flag_isLoaded = 0;
    end
    
    
    %----------------------------------------------------------------------
    % ---------------------------------------------------------------------
    %   End: Loop over load step
    % ---------------------------------------------------------------------
    %----------------------------------------------------------------------
    
    
    fprintf('\n');
    fprintf('  End of the problem.\n\n');
    fprintf('----------------------------------------------------------------\n');
    fprintf('----------------------------------------------------------------\n\n');
end
