%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  This code is based on that by Mike Borden in the Bezier extraction paper.
%  
%  
%  Summary:
%  
%  This routine computes the element Bezier extraction matrices for a
%  given knot vector in 1D (this is assumed to be an open knot vector).
%  We use these matrices to map the Bernstein polynomial basis functions
%  to the B-spline basis functions restricted to an element.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open.
%  
%  
%  Instructions:
%  
%  Type one of the following onto Matlab's command window or in a code,
%  
%      [bezierExtractions, nodeIndexShifts, numElements, elementSizes] = build_bezier_extraction(knots, p);
%      [bezierExtractions, nodeIndexShifts, numElements] = build_bezier_extraction(knots, p);
%      [bezierExtractions, nodeIndexShifts] = build_bezier_extraction(knots, p);
%      bezierExtractions = build_bezier_extraction(knots, p);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      p is the degree of the B-splines
%  
%  
%  Outputs:
%  
%  1. (p + 1) x (p + 1) x (numElements) array of Bezier extraction matrices
%  
%      bezierExtractions(:, :, e) returns the Bezier extraction matrix for
%      the e-th element
%  
%  2. (numElements) x 1 column vector of shifts
%  
%      The shifts are added to the node indices e : (e + p), which is
%      applicable to a knot vector with no repeated interior knots.
%  
%  3. Number of elements
%  
%  4. (numElements) x 1 column vector of element sizes
%  
%      The element sizes are used as Jacobians to get the derivatives of
%      the B-splines from the derivatives of the Bernstein polynomials.
%--------------------------------------------------------------------------
function [bezierExtractions, nodeIndexShifts, numElements, elementSizes] = build_bezier_extraction(knots, p)
    % Some useful constants
    constant_pp1 = p + 1;
    
    numNodes = size(knots, 1) - constant_pp1;
    numElements = size(unique(knots), 1) - 1;
    
    
    %----------------------------------------------------------------------
    %  Initialize the extraction matrices, the node index shift vector,
    %  and the element sizes vector
    %----------------------------------------------------------------------
    % Initialize the extraction matrices to be the identity matrix
    bezierExtractions = zeros(constant_pp1, constant_pp1, numElements);
    
    for e = 1 : numElements
        bezierExtractions(:, :, e) = eye(constant_pp1);
    end
    
    % Initially, we assume that none of the interior knots are repeated
    nodeIndexShifts = zeros(numElements, 1);
    
    % Element sizes
    elementSizes = zeros(numElements, 1);
    
    
    %----------------------------------------------------------------------
    %  Initialize variables for the Bezier element [xi1, xi2]
    %----------------------------------------------------------------------
    % Set xi1 to the first knot
    xi1 = knots(1);
    
    % Initialize the coefficients for knot insertion
    alphas_knotInsertion = zeros(p - 1, 1);
    
    % Indices of the knots corresponding to the first Bezier element
    a = constant_pp1;
    b = a + 1;
    
    % Element index
    e = 1;
    
    % Keep track of the total number of shifts
    shift_total = 0;
    
    
    %----------------------------------------------------------------------
    %  Loop over the knots
    %----------------------------------------------------------------------
    while (b <= numNodes)
        % Count how many times the knot xi2 is repeated
        i = b;
        
        while (b <= numNodes && knots(b) == knots(b + 1))
            b = b + 1;
        end
        
        xi2 = knots(b);
        
        numRepeats = b - i + 1;
        
        % Update the shift if there are repeated knots
        shift_total = shift_total + (numRepeats - 1);
        
        % Record the element size
        elementSizes(e) = xi2 - xi1;
        
        
        %------------------------------------------------------------------
        %  Knot insertion
        %------------------------------------------------------------------
        % Insert the knot xi2 to get a Bezier element
        numKnotsToInsert = p - numRepeats;
        
        if (numKnotsToInsert > 0)
            % Compute the coefficients for the knot insertion
            numerator = xi2 - xi1;
            
            for i = 1 : numKnotsToInsert
                alphas_knotInsertion(i) = numerator / (knots(a + numRepeats + i) - xi1);
            end
            
            % Compute the nodes for the Bezier element
            for i = 1 : numKnotsToInsert
                s = numRepeats + i;
                
                for k = constant_pp1 : -1 : (s + 1)
                    alpha = alphas_knotInsertion(k - s);
                    
                    bezierExtractions(:, k, e) = alpha * bezierExtractions(:, k, e) + (1 - alpha) * bezierExtractions(:, k - 1, e);
                end
                
                % Update the overlapping entries of the next extraction matrix
                if (b <= numNodes)
                    index = numKnotsToInsert - i + 1;
                    bezierExtractions((index : (index + i))', index, e + 1) = bezierExtractions(((constant_pp1 - i) : constant_pp1)', constant_pp1, e);
                end
            end
            
            % Update indices for the next extraction matrix
            if (b <= numNodes)
                a = b;
                b = a + 1;
                
                xi1 = xi2;
            end
            
            e = e + 1;
            nodeIndexShifts(e) = shift_total;
            
        elseif (numKnotsToInsert == 0)
            % Update indices for the next extraction matrix
            if (b <= numNodes)
                a = b;
                b = a + 1;
                
                xi1 = xi2;
                
                e = e + 1;
                nodeIndexShifts(e) = shift_total;
            end
        end
    end
    
    % Record the element size for the last element
    elementSizes(e) = knots(b) - xi1;
    
    % Record the number of elements for future use
    numElements = e;
end