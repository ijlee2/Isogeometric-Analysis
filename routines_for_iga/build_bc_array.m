%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  
%  Summary:
%  
%  This routine builds the BCU (displacement BC) and BCF (force BC) arrays
%  from the displacement and force BC arrays.
%  
%  
%  Warning:
%  
%  Use this routine after calling build_id_array routine and creating the
%  BCs_displacement and BCs_force arrays. These arrays have 3 columns of
%  the form,
%  
%      [global node index, dof index, BC value]
%  
%  Note that we will assume zero force for all DOFs on which displacement
%  or force BC has not been specified. This allows us to specify only the 
%  nonzero forces in the BCs_force array and save memory.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      [BCU_array, BCF_array] = build_bc_array(BCs_displacement, BCs_force, ID_array);
%  
%  where,
%  
%      BCs_displacement, BCs_force are the BC arrays
%      ID_array is the ID ("Destination") array
%  
%  
%  Output:
%  
%  1. (numBCs_displacement) x 2 array
%  
%      Prescribe the displacement values by u(BCU_array(:, 1)) = BCU_array(:, 2);
%  
%  2. (numBCs_force) x 2 array
%  
%      Prescribe the force values by f(BCF_array(:, 1)) = BCF_array(:, 2);
%--------------------------------------------------------------------------
function [BCU_array, BCF_array] = build_bc_array(BCs_displacement, BCs_force, ID_array)
    %----------------------------------------------------------------------
    %  Remove duplicate BC information
    %----------------------------------------------------------------------
    BCs_displacement = unique(BCs_displacement, 'rows');
    BCs_force = unique(BCs_force, 'rows');
    
    numBCs_displacement = size(BCs_displacement, 1);
    numBCs_force = size(BCs_force, 1);
    
    
    % Check for errors
    if (numBCs_displacement > 0 && numBCs_displacement ~= size(unique(BCs_displacement(:, [1 2]), 'rows'), 1))
        fprintf('Warning: The displacement BCs array contains two identical DOFs with different values.\n\nNo operation has been performed.\n\n');
        
        BCU_array = [];
        BCF_array = [];
        
        return;
    end
    
    if (numBCs_force > 0 && numBCs_force ~= size(unique(BCs_force(:, [1 2]), 'rows'), 1))
        fprintf('Warning: The force BCs array contains two identical DOFs with different values.\n\nNo operation has been performed.\n\n');
        
        BCU_array = [];
        BCF_array = [];
        
        return;
    end
    
    
    %----------------------------------------------------------------------
    %  Create the BCU array
    %----------------------------------------------------------------------
    BCU_array = zeros(numBCs_displacement, 2);
    
    for i = 1 : numBCs_displacement
        BCU_array(i, 1) = ID_array(BCs_displacement(i, 2), BCs_displacement(i, 1));
        BCU_array(i, 2) = BCs_displacement(i, 3);
    end
    
    
    %----------------------------------------------------------------------
    %  Create the BCF array
    %----------------------------------------------------------------------
    BCF_array = zeros(numBCs_force, 2);
    
    for i = 1 : numBCs_force
        BCF_array(i, 1) = ID_array(BCs_force(i, 2), BCs_force(i, 1));
        BCF_array(i, 2) = BCs_force(i, 3);
    end
end