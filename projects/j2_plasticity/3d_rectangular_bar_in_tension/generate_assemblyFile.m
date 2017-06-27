%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine creates the assembly files for the problem of a 3D
%  rectangular bar in tension. We use quadratic NURBS as basis functions.
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
%      numRefinements is the number of times we h-refine the knot vectors
%  
%  
%  Output:
%  
%  1. Assembly files (.mat files)
%--------------------------------------------------------------------------
function generate_assemblyFile(numRefinements)
    % Set the path to the assembly files directory
    directory_assembly = sprintf('../assembly_files/numRefinements%d/', numRefinements);
    
    % Create the directory if it does not exist
    if ~exist(directory_assembly, 'dir')
        mkdir(directory_assembly);
    end
    
    % File paths
    file_assembly_global = strcat(directory_assembly, 'file_assembly_global');
    file_assembly_patch  = strcat(directory_assembly, 'file_assembly_patch');
    
    
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
    numDOFsPerNode = 3;
    
    % Shared nodes array
    SN_array = [];
    
    % Global nodes array
    GN_array = [];
    
    
    %----------------------------------------------------------------------
    %  Initialize parameters for plasticity
    %----------------------------------------------------------------------
    % Yield stress (Pa)
    material_sigmaY = 500e6;
    
    % Hardening modulus (Pa)
    material_H = 550e6;
    
    % Isotropic hardening factor (between 0 and 1)
    isotropicHardeningFactor = 1;
    
    % Number of load steps
    numLoadSteps = 1;
    
    % Number of points where stresses and strains are evaluated
    numStressPoints = 0;
    
    % Frequency of saves
    numTimeStepsBetweenSaves = 1;
    
    % Parameters for Newton's method
    maxNewtonsMethod = 10;
    tolNewtonsMethod = 1e-14;
    
    
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
        struct('elementType'        , 'NURBS', ...
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
    patches(1).description = strcat('3D rectangular bar specimen, refinement level', sprintf(' %d', numRefinements));
    
    
    %----------------------------------------------------------------------
    %  Set the material parameters
    %----------------------------------------------------------------------
    patches(1).material_E = 210e9;    % Young's modulus (Pa)
    patches(1).material_nu = 0.3;     % Poisson's ratio (dimensionless)
    
    
    %----------------------------------------------------------------------
    %  Set the degrees
    %----------------------------------------------------------------------
    patches(1).p1 = 2;
    patches(1).p2 = 1;
    patches(1).p3 = 1;
    
    
    %----------------------------------------------------------------------
    %  Set the knots
    %----------------------------------------------------------------------
    patches(1).knots1 = [0; 0; 0; 1; 1; 2; 3; 4; 5; 5; 6; 6; 6];
    patches(1).knots2 = [0; 0; 1; 1];
    patches(1).knots3 = [0; 0; 1; 1];
    
    
    %----------------------------------------------------------------------
    %  Set the nodes
    %----------------------------------------------------------------------
    % Parameters for specimen size (all lengths in m)
    grip_length = 75e-3;              % Length of the grip
    grip_width = 50e-3;               % Width of the grip
    
    fillet_length = 15e-3;            % Length of the fillet
    fillet_width = 5e-3;              % Width of the fillet
    fillet_radius = sqrt(((fillet_length^2 - fillet_width^2) / (2 * fillet_width))^2 + fillet_length^2);
    w_fillet = sqrt(1 - fillet_width / (2 * fillet_radius));
    
    bar_length = 220e-3;              % Length of the bar
%   bar_width = 40e-3;                % Width of the bar
%   gauge_length = 200e-3;            % Length of the gauge section
    
    bar_totalLength = grip_length + fillet_length + bar_length + fillet_length + grip_length;
    bar_thickness = 5e-3;             % Thickness of the bar (m)
    
    x_loc = [0; ...
             grip_length / 2; ...
             grip_length; ...
             grip_length; ...
             grip_length + fillet_length; ...
             grip_length + fillet_length + bar_length; ...
             grip_length + fillet_length + bar_length + fillet_length; ...
             grip_length + fillet_length + bar_length + fillet_length; ...
             grip_length + fillet_length + bar_length + fillet_length + grip_length / 2; ...
             bar_totalLength];
    
    y_loc = [0; ...
             0; ...
             0; ...
             fillet_width; ...
             fillet_width; ...
             fillet_width; ...
             fillet_width; ...
             0; ...
             0; ...
             0];
    
    z_loc = [0; ...
             0; ...
             0; ...
             0; ...
             0; ...
             0; ...
             0; ...
             0; ...
             0; ...
             0];
    
    w_loc = [1; ...
             1; ...
             1; ...
             w_fillet; ...
             1; ...
             1; ...
             w_fillet; ...
             1; ...
             1; ...
             1];
    
    patches(1).nodes = [x_loc, y_loc             , z_loc                , w_loc; ...
                        x_loc, grip_width - y_loc, z_loc                , w_loc; ...
                        x_loc, y_loc             , bar_thickness - z_loc, w_loc; ...
                        x_loc, grip_width - y_loc, bar_thickness - z_loc, w_loc];
    
    
    %----------------------------------------------------------------------
    %  Project the NURBS nodes to B-spline nodes for degree elevation and
    %  knot refinement
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        patches(p).nodes = project_up(patches(p).nodes);
    end
    
    
    %----------------------------------------------------------------------
    %  Elevate the degrees
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        % Elevate the degree to 2 in all directions
        [patches(p).knots1, patches(p).knots2, patches(p).knots3, patches(p).nodes, patches(p).p1, patches(p).p2, patches(p).p3] = ...
            refine_p_volume(patches(p).knots1, ...
                            patches(p).knots2, ...
                            patches(p).knots3, ...
                            patches(p).nodes, ...
                            patches(p).p1, ...
                            patches(p).p2, ...
                            patches(p).p3, ...
                            0, 1, 1);
    end
    
    
    %----------------------------------------------------------------------
    %  Refine the knots
    %----------------------------------------------------------------------
    % Specify the knots for insertion
    switch numRefinements
        case 0
            temp = 0.5 + 0.5 * sqrt(setdiff((0 : 0.5^3 : 1)', [0; 1]));
            
            knotsForInsertion1 = [0.5; 2.45; 2.6; (2 + temp); sort(4 - temp); 3.4; 3.55; 5.5];
            knotsForInsertion2 = [0.25; 0.5; 0.75];
            knotsForInsertion3 = [0.25; 0.5; 0.75];
            
        case 1
            temp = 0.5 + 0.5 * sqrt(setdiff((0 : 0.5^4 : 1)', [0; 1]));
            
            knotsForInsertion1 = [0.25; 0.5; 0.75; 2.3; 2.45; 2.55; (2 + temp); sort(4 - temp); 3.45; 3.55; 3.7; 5.25; 5.5; 5.75];
            knotsForInsertion2 = [0.25; 0.5; 0.75];
            knotsForInsertion3 = [0.25; 0.5; 0.75];
            
        case 2
            temp = 0.5 + 0.5 * sqrt(setdiff((0 : 0.5^5 : 1)', [0; 1]));
            
            knotsForInsertion1 = [0.25; 0.5; 0.75; 2.3; 2.45; 2.55; (2 + temp); sort(4 - temp); 3.45; 3.55; 3.7; 5.25; 5.5; 5.75];
            knotsForInsertion2 = [0.125; 0.25; 0.375; 0.5; 0.625; 0.75; 0.875];
            knotsForInsertion3 = [0.25; 0.5; 0.75];
            
        otherwise
            fprintf('Error: Please specify the knots for insertion.\n\n');
            
            return;
            
    end
    
    clear temp
    
    
    % Refine the knot vectors
    [patches(1).knots1, patches(1).knots2, patches(1).knots3, patches(1).nodes] = ...
        refine_h_volume(patches(1).knots1, ...
                        patches(1).knots2, ...
                        patches(1).knots3, ...
                        patches(1).nodes, ...
                        patches(1).p1, ...
                        patches(1).p2, ...
                        patches(1).p3, ...
                        knotsForInsertion1, knotsForInsertion2, knotsForInsertion3);
    
    clear knotsForInsertion1 knotsForInsertion2 knotsForInsertion3;
    
    
    %----------------------------------------------------------------------
    %  Project the B-spline nodes back to NURBS nodes
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        patches(p).nodes = project_down(patches(p).nodes);
        
        %{
        % Draw the NURBS patch (for debugging)
        draw_nurbs_volume(patches(p).knots1, patches(p).knots2, patches(p).knots3, patches(p).nodes, patches(p).p1, patches(p).p2, patches(p).p3);
        axis image;
        return;
        %}
    end
    
    
    %----------------------------------------------------------------------
    %  Build the Bezier extractions and the IEN array
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        [patches(p).bezierExtractions1, nodeIndexShifts1, patches(p).numElements1, patches(p).elementSizes1] = build_bezier_extraction(patches(p).knots1, patches(p).p1);
        [patches(p).bezierExtractions2, nodeIndexShifts2, patches(p).numElements2, patches(p).elementSizes2] = build_bezier_extraction(patches(p).knots2, patches(p).p2);
        [patches(p).bezierExtractions3, nodeIndexShifts3, patches(p).numElements3, patches(p).elementSizes3] = build_bezier_extraction(patches(p).knots3, patches(p).p3);
        
        patches(p).IEN_array = build_ien_array(patches(p).knots1, ...
                                               patches(p).knots2, ...
                                               patches(p).knots3, ...
                                               patches(p).p1, ...
                                               patches(p).p2, ...
                                               patches(p).p3, ...
                                               nodeIndexShifts1, ...
                                               nodeIndexShifts2, ...
                                               nodeIndexShifts3);
    end
    
    clear nodeIndexShifts1 nodeIndexShifts2 nodeIndexShifts3;
    
    
    %----------------------------------------------------------------------
    %  Find the numbers of nodes, elements, and matrix entries
    %----------------------------------------------------------------------
    for p = 1 : numPatches
        % Number of nodes along each direction
        patches(p).numNodes1 = size(patches(p).knots1, 1) - (patches(p).p1 + 1);
        patches(p).numNodes2 = size(patches(p).knots2, 1) - (patches(p).p2 + 1);
        patches(p).numNodes3 = size(patches(p).knots3, 1) - (patches(p).p3 + 1);
        
        % Number of nodes on the patch
        numNodesOnPatch(p) = patches(p).numNodes1 * patches(p).numNodes2 * patches(p).numNodes3;
        
        % Number of nodes before the patch
        if (p > 1)
            numNodesBeforePatch(p) = numNodesBeforePatch(p - 1) + numNodesOnPatch(p - 1);
        end
        
        % Number of elements on the patch
        numElementsOnPatch(p) = patches(p).numElements1 * patches(p).numElements2 * patches(p).numElements3;
        
        % Number of nodes on each element along each direction
        patches(p).numNodesPerElement1 = patches(p).p1 + 1;
        patches(p).numNodesPerElement2 = patches(p).p2 + 1;
        patches(p).numNodesPerElement3 = patches(p).p3 + 1;
        patches(p).numNodesPerElement = patches(p).numNodesPerElement1 * patches(p).numNodesPerElement2 * patches(p).numNodesPerElement3;
        
        % Number of DOFs on each element
        patches(p).numDOFsPerElement = numDOFsPerNode * patches(p).numNodesPerElement;
        
        % Number of matrix entries that we will compute as a result
        numMatrixEntriesPerElement = patches(p).numDOFsPerElement^2;
        numMatrixEntries = numMatrixEntries + numMatrixEntriesPerElement * numElementsOnPatch(p);
    end
    
    
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
    for p = 1 : numPatches
        patches(p).numQuadraturePoints = [3; 3; 3];
        
        numStressPoints = numStressPoints + prod(patches(p).numQuadraturePoints) * numElementsOnPatch(p);
    end
    
    
    %----------------------------------------------------------------------
    %  Save the assembly files
    %----------------------------------------------------------------------
    % Assembly file for the global structure
    save(file_assembly_global, ...
         'numPatches'              , ...
         'numElementsOnPatch'      , ...
         'numNodesOnPatch'         , ...
         'numNodesBeforePatch'     , ...
         'numMatrixEntries'        , ...
         'numDOFs'                 , ...
         'numDOFsPerNode'          , ...
         'SN_array'                , ...
         'GN_array'                , ...
         ...
         'material_sigmaY'         , ...
         'material_H'              , ...
         'isotropicHardeningFactor', ...
         'numLoadSteps'            , ...
         'numStressPoints'         , ...
         'numTimeStepsBetweenSaves', ...
         'maxNewtonsMethod'        , ...
         'tolNewtonsMethod'        , ...
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
    generate_bcFile(numRefinements);
end
