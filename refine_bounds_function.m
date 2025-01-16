function [E, R] = refine_bounds_function( A, b)
% CONSTRAINT_BOUNDS - Compute bounds for constraints in a zonotope representation.
%
% Inputs:

%   A - Constraint matrix (nc x ng)
%   b - Constraint vector (nc x 1)

%
% Outputs:
%   E - Final bounds matrix for generators (ng x 2)
%   R - Final bounds matrix for constraints (ng x 2)


    % Dimensions
    [nc, ng] = size(A); % nc = number of constraints, ng = number of generators

    % Initialize bounds
    E = repmat([-1, 1], ng, 1); % E initialized as [-1, 1] for all generators
    R = repmat([-inf, inf], ng, 1); % R initialized as [-Inf, Inf] for all generators

        % Iterate over constraints and generators
        for i = 1:nc
            for j = 1:ng
                ajj = A(i, j); % Coefficient of xi_j in constraint i
                if ajj == 0
                    continue; % Skip if xi_j is not in this constraint
                end

                % Calculate bounds from constraint i
                sum_rest_low = 0; % Contribution of other variables (k != j)
                sum_rest_high = 0; % Contribution of other variables (k != j)

                for k = 1:ng
                    if k ~= j
                        R_temp_low = A(i, k) / ajj * E(k, 1);
                        R_temp_high = A(i, k) / ajj * E(k, 2);

                        if R_temp_low > R_temp_high
                            temp = R_temp_low;
                            R_temp_low = R_temp_high;
                            R_temp_high = temp;
                        end

                        sum_rest_low = sum_rest_low + R_temp_low; % Use lower bound of Ek                       
                        sum_rest_high = sum_rest_high + R_temp_high; % Use upper bound of Ek

                    end
                end

                sum_rest = [sum_rest_low, sum_rest_high];
                Rj = b(i) / ajj - sum_rest;           

                % Compute Rj bounds
                Rj_lower = Rj(1);
                Rj_upper = Rj(2);

                if Rj_lower > Rj_upper
                    a = Rj_lower;
                    Rj_lower = Rj_upper;
                    Rj_upper = a;
                end

                % Update Rj
                R(j, 1) = max(R(j, 1), Rj_lower);
                R(j, 2) = min(R(j, 2), Rj_upper);

                E(j, 1) = max(E(j, 1), R(j, 1)); % Update lower bound
                E(j, 2) = min(E(j, 2), R(j, 2)); % Update upper bound
                % 
                if E(j, 1) > E(j, 2)
                    a = E(j, 1);
                    E(j, 1) = E(j, 2);
                    E(j, 2) = a;
                end
            end
        end
    end

