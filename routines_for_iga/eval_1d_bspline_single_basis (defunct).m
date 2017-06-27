%--------------------------------------------------------------------------
%  Author: Isaac J. Lee (ijlee2@ices.utexas.edu)
%  
%  This code is based on that by Piegl and Tiller in their NURBS book.
%  
%  
%  Summary:
%  
%  This routine evaluates the 1D B-spline of degree p at multiple locations.
%  
%  
%  Warning:
%  
%  The knot vector is assumed to be open.
%  
%  This routine will be inefficient for finite element analysis. Please
%  use bezier_extraction and eval_bernstein_all to compute the Bernstein
%  polynomials, and find the B-splines using the affine map between the
%  the Bernstein and B-spline domains.
%  
%  
%  Instructions:
%  
%  Type the following onto Matlab's command window or in a code,
%  
%      output = eval_1d_bspline_single_basis(knots, xi, i, p);
%  
%  where,
%  
%      knots is the knot vector (column vector)
%      xi is an array of knot values (column vector)
%      i is the index of the B-spline
%      p is the degree of the B-spline
%  
%  
%  Output:
%  
%  1. (output_size) x 1 column vector
%  
%      output returns the values of the B-spline at the knots.
%--------------------------------------------------------------------------
function output = eval_1d_bspline_single_basis(knots, xi, i, p)
    output_size = size(xi, 1);
    output = zeros(output_size, 1);
    
    for index = 1 : output_size
        % If the knot is open, then we automatically set the evaluation at
        % the last knot value for the last basis function to be 1
        if (i == output_size - p - 1 && xi(index) == knots(output_size))
            output(index) = 1;
            
        % If xi lies outside the support of the i-th basis function, then
        % we automatically set the evaluation to be 0
        elseif (xi(index) < knots(i) || xi(index) >= knots(i + p + 1))
            % Nothing to be done
            
        % Evaluate N_{i, p} at xi using the Cox-de Boor recursion
        else
            N = zeros(p + 1, 1);
            
            % Evaluate the B-spline basis functions N_{i+j, 0} at xi
            for j = 0 : p
                if (knots(i + j) <= xi(index) && xi(index) < knots(i + j + 1))
                    N(j + 1) = 1;
                else
                    N(j + 1) = 0;
                end
            end
            
            % Compute the triangular table
            for k = 1 : p
                if (N(1) == 0)
                    saved = 0;
                else
                    saved = (xi(index) - knots(i))/(knots(i + k) - knots(i)) * N(1);
                end
                
                for j = 1 : p - k + 1
                    Xi_left = knots(i + j);
                    Xi_right = knots(i + j + k);
                    
                    if (N(j + 1) == 0)
                        N(j) = saved;
                        saved = 0;
                    else
                        temp = N(j + 1)/(Xi_right - Xi_left);
                        N(j) = saved + (Xi_right - xi(index))*temp;
                        saved = (xi(index) - Xi_left)*temp;
                    end
                end
            end
            
            output(index) = N(1);
        end
    end
end