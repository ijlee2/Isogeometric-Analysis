%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine creates the assembly file for the problem of 1D cracked bar
%  in tension. We use B-splines as basis functions.
%  
%  
%  Instructions:
%  
%  Use the terminal to run the driver file with this command:
%  
%      ./matbg.sh driver_generate_assemblyFile.m output
%  
%  Note that,
%  
%      order is the order of the phase field theory (2, 4, 6, 8)
%      numRefinements is the number of times we h-refine the knot vectors
%  
%  
%  Output:
%  
%  1. Assembly files (.mat files)
%--------------------------------------------------------------------------
function generate_assemblyFile(order, numRefinements)
    % Check for errors
    if (order ~= 2 && order ~= 4 && order ~= 6 && order ~= 8)
        fprintf('Error: The order of the phase field theory must be 2, 4, 6, or 8.\n\n');
        
        return;
    end
    
    
    % Set the path to the assembly files directory
    path_to_assembly_directory = sprintf('../assembly_files/order%d/numRefinements%d/', order, numRefinements);
    
    % Create the directory if it does not exist
    if ~exist(path_to_assembly_directory, 'dir')
        mkdir(path_to_assembly_directory);
    end
    
    % File paths
    file_assembly_global = strcat(path_to_assembly_directory, 'file_assembly_global');
    file_assembly_patch  = strcat(path_to_assembly_directory, 'file_assembly_patch');
    
    
    %----------------------------------------------------------------------
    %  Initialize parameters for FEM
    %----------------------------------------------------------------------
    % Number of patches
    numPatches = 1;
    
    % Number of elements on each patch
    numElementsOnPatch = zeros(numPatches, 1);
    
    % Number of nodes on each patch
    numNodesOnPatch = zeros(numPatches, 1);
    
    % Number of nodes encountered before the patch
    numNodesBeforePatch = zeros(numPatches, 1);
    
    % Number of matrix entries that will be computed
    numMatrixEntries = 0;
    
    % Number of degrees of freedom (DOFs)
    numDOFs = 0;
    
    % Number of DOFs (fields) defined on each node
    numDOFsPerNode = 1;
    
    % Shared nodes array
    SN_array = [];
    
    % Global nodes array
    GN_array = [];
    
    
    %----------------------------------------------------------------------
    %  Initialize parameters for 1D phase field theory
    %----------------------------------------------------------------------
    % Half-length of the bar (m)
    material_L = 100;
    
    % Cross-sectional area (m^2)
    material_A = 1;
    
    % Critical energy release rate, per unit area (J/m^2 = N/m)
    material_G_c = 1;
    
    % Phase field parameter (m)
    material_ell_0 = 1/8;
    
    % Maximum number of alternations between solving the momentum and
    % phase field equations
    maxAlternations = 10;
    
    
    %----------------------------------------------------------------------
    %  Create an array of NURBS structs
    %  
    %  description
    %      -- the description of the NURBS patch used for debugging 
    %      -- string
    %  
    %  knots1, knots2, knots3
    %      -- the knot vectors in directions 1, 2, and 3
    %      -- column vectors
    %  
    %  nodes
    %      -- the control points in the physical space
    %      -- matrix with d columns, where d is the dimension
    %  
    %  p1, p2, p3
    %      -- the degree of the NURBS in directions 1, 2, and 3
    %      -- positive integers
    %----------------------------------------------------------------------
    patches(numPatches) = ...
        struct('elementType'        , 'B-spline', ...
               'description'        , '', ...
               'material_E'         , 0, ...
               'material_nu'        , 0, ...
               'p1'                 , 0, ...
               'p2'                 , 0, ...
               'p3'                 , 0, ...
               'knots1'             , [], ...
               'knots2'             , [], ...
               'knots3'             , [], ...
               'nodes'              , [], ...
               'numNodes1'          , 0, ...
               'numNodes2'          , 0, ...
               'numNodes3'          , 0, ...
               'numNodesPerElement1', 0, ...
               'numNodesPerElement2', 0, ...
               'numNodesPerElement3', 0, ...
               'numNodesPerElement' , 0, ...
               'numDOFsPerElement'  , 0, ...
               'bezierExtractions1' , [], ...
               'bezierExtractions2' , [], ...
               'bezierExtractions3' , [], ...
               'numElements1'       , 0, ...
               'numElements2'       , 0, ...
               'numElements3'       , 0, ...
               'elementSizes1'      , [], ...
               'elementSizes2'      , [], ...
               'elementSizes3'      , [], ...
               'IEN_array'          , [], ...
               'numQuadraturePoints', []);
    
    
    %----------------------------------------------------------------------
    %  Set the description
    %----------------------------------------------------------------------
    patches(1).description = strcat('1D bar with a crack in the middle, refinement level', sprintf(' %d', numRefinements));
    
    
    %----------------------------------------------------------------------
    %  Set the material parameters
    %----------------------------------------------------------------------
    % Young's modulus (Pa)
    patches(1).material_E = 1e6;
    
    % Poisson's ratio (dimensionless)
    patches(1).material_nu = 0.3;
    
    
    %----------------------------------------------------------------------
    %  Set the degrees
    %----------------------------------------------------------------------
    patches(1).p1 = 1;
    
    
    %----------------------------------------------------------------------
    %  Set the knots
    %----------------------------------------------------------------------
    patches(1).knots1 = [0; 0; 1; 1];
    
    
    %----------------------------------------------------------------------
    %  Set the nodes
    %----------------------------------------------------------------------
    patches(1).nodes = [-material_L; ...
                         material_L];
    
    
    %----------------------------------------------------------------------
    %  Project the NURBS nodes to B-spline nodes for degree elevation and
    %  knot refinement
    %----------------------------------------------------------------------
    % Do nothing
    
    
    %----------------------------------------------------------------------
    %  Elevate the degrees
    %----------------------------------------------------------------------
    % Elevate the degree to (order / 2) in the x-direction
    if (patches(1).p1 < (order / 2))
        [patches(1).knots1, patches(1).nodes, patches(1).p1] = ...
            refine_p_curve(patches(1).knots1, ...
                           patches(1).nodes, ...
                           patches(1).p1, ...
                           (order / 2) - patches(1).p1);
    end
    
    
    %----------------------------------------------------------------------
    %  Refine the knots
    %----------------------------------------------------------------------
    flag_uniformlyRefine = 0;
    
    % Uniform refinement
    if (flag_uniformlyRefine == 1)
        % Set the desired number of elements
        numElements = 4000 * 2^numRefinements;
        
        % Set the knots for insertion
        knotsForInsertion1 = setdiff((0 : 1 / numElements : 1)', [0; 1]);
        
    % Non-uniform refinement
    else
        % Set the desired number of elements
        numElements = 1000 * 2^numRefinements;
        
        % Set the knots for insertion
        temp = 0.5 * (setdiff((0 : 1 / numElements : 1)', [0; 1])).^(1/4);
        
        knotsForInsertion1 = [temp; 0.5; sort(1 - temp)];
        
        clear temp;
        
    end
    
    % Refine the knot vectors
    [patches(1).knots1, patches(1).nodes] = ...
        refine_h_curve(patches(1).knots1, ...
                       patches(1).nodes, ...
                       patches(1).p1, ...
                       knotsForInsertion1);
    
    clear knotsForInsertion1;
    
    
    %----------------------------------------------------------------------
    %  Project the B-spline nodes back to NURBS nodes
    %----------------------------------------------------------------------
    % Do nothing
    
    
    %----------------------------------------------------------------------
    %  Build the Bezier extractions and the IEN array
    %----------------------------------------------------------------------
    [patches(1).bezierExtractions1, nodeIndexShifts1, patches(1).numElements1, patches(1).elementSizes1] = build_bezier_extraction(patches(1).knots1, patches(1).p1);
    
    patches(1).IEN_array = build_ien_array(patches(1).knots1, ...
                                           [], ...
                                           [], ...
                                           patches(1).p1, ...
                                           [], ...
                                           [], ...
                                           nodeIndexShifts1, ...
                                           [], ...
                                           []);
    
    clear nodeIndexShifts1;
    
    
    %----------------------------------------------------------------------
    %  Find the numbers of nodes, elements, and matrix entries
    %----------------------------------------------------------------------
    % Number of nodes along each direction
    patches(1).numNodes1 = size(patches(1).knots1, 1) - (patches(1).p1 + 1);
    
    % Number of nodes on the patch
    numNodesOnPatch(1) = patches(1).numNodes1;
    
    % Number of nodes before the patch (do nothing, since we have 1 patch)
    
    % Number of elements on the patch
    numElementsOnPatch(1) = patches(1).numElements1;
    
    % Number of nodes on each element along each direction
    patches(1).numNodesPerElement1 = patches(1).p1 + 1;
    patches(1).numNodesPerElement = patches(1).numNodesPerElement1;
    
    % Number of DOFs on each element
    patches(1).numDOFsPerElement = numDOFsPerNode * patches(1).numNodesPerElement;
    
    % Number of matrix entries that we will compute as a result
    numMatrixEntriesPerElement = patches(1).numDOFsPerElement^2;
    numMatrixEntries = numMatrixEntries + numMatrixEntriesPerElement * numElementsOnPatch(1);
    
    
    %----------------------------------------------------------------------
    %  Build the SN array
    %  [my patch index, my node index, target patch index, target node index]
    %----------------------------------------------------------------------
    SN_array = [];
    
    
    %----------------------------------------------------------------------
    %  Build the GN array
    %----------------------------------------------------------------------
    GN_array = build_gn_array(SN_array, numNodesOnPatch);
    
    % Number of degrees of freedom
    numDOFs = numDOFsPerNode * size(unique(GN_array), 1);
    
    
    %----------------------------------------------------------------------
    %  Set the quadrature rule
    %----------------------------------------------------------------------
    patches(1).numQuadraturePoints = order / 2 + 1;
    
    
    %----------------------------------------------------------------------
    %  Save the assembly files
    %----------------------------------------------------------------------
    % Assembly file for the global structure
    save(file_assembly_global, ...
         'numPatches'         , ...
         'numElementsOnPatch' , ...
         'numNodesOnPatch'    , ...
         'numNodesBeforePatch', ...
         'numMatrixEntries'   , ...
         'numDOFs'            , ...
         'numDOFsPerNode'     , ...
         'SN_array'           , ...
         'GN_array'           , ...
         ...
         'material_L'         , ...
         'material_A'         , ...
         'material_G_c'       , ...
         'material_ell_0'     , ...
         ...
         'maxAlternations'    , ...
         '-v7.3');
    
    % Assembly file for the patches
    for p = 1 : numPatches
        elementType         = patches(p).elementType;
        description         = patches(p).description;
        material_E          = patches(p).material_E;
        material_nu         = patches(p).material_nu;
        p1                  = patches(p).p1;
        p2                  = patches(p).p2;
        p3                  = patches(p).p3;
        knots1              = patches(p).knots1;
        knots2              = patches(p).knots2;
        knots3              = patches(p).knots3;
        nodes               = patches(p).nodes;
        numNodes1           = patches(p).numNodes1;
        numNodes2           = patches(p).numNodes2;
        numNodes3           = patches(p).numNodes3;
        numNodesPerElement1 = patches(p).numNodesPerElement1;
        numNodesPerElement2 = patches(p).numNodesPerElement2;
        numNodesPerElement3 = patches(p).numNodesPerElement3;
        numNodesPerElement  = patches(p).numNodesPerElement;
        numDOFsPerElement   = patches(p).numDOFsPerElement;
        bezierExtractions1  = patches(p).bezierExtractions1;
        bezierExtractions2  = patches(p).bezierExtractions2;
        bezierExtractions3  = patches(p).bezierExtractions3;
        numElements1        = patches(p).numElements1;
        numElements2        = patches(p).numElements2;
        numElements3        = patches(p).numElements3;
        elementSizes1       = patches(p).elementSizes1;
        elementSizes2       = patches(p).elementSizes2;
        elementSizes3       = patches(p).elementSizes3;
        IEN_array           = patches(p).IEN_array;
        numQuadraturePoints = patches(p).numQuadraturePoints;
        
        save(sprintf('%s%d', file_assembly_patch, p), ...
             'elementType'        , ...
             'description'        , ...
             'material_E'         , ...
             'material_nu'        , ...
             'p1'                 , ...
             'p2'                 , ...
             'p3'                 , ...
             'knots1'             , ...
             'knots2'             , ...
             'knots3'             , ...
             'nodes'              , ...
             'numNodes1'          , ...
             'numNodes2'          , ...
             'numNodes3'          , ...
             'numNodesPerElement1', ...
             'numNodesPerElement2', ...
             'numNodesPerElement3', ...
             'numNodesPerElement' , ...
             'numDOFsPerElement'  , ...
             'bezierExtractions1' , ...
             'bezierExtractions2' , ...
             'bezierExtractions3' , ...
             'numElements1'       , ...
             'numElements2'       , ...
             'numElements3'       , ...
             'elementSizes1'      , ...
             'elementSizes2'      , ...
             'elementSizes3'      , ...
             'IEN_array'          , ...
             'numQuadraturePoints', ...
             '-v7.3');
    end
    
    
    %----------------------------------------------------------------------
    %  Set the boundary conditions
    %----------------------------------------------------------------------
    generate_bcFile(order, numRefinements);
end
