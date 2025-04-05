% load conzono_variable.mat
load conzono_MoreComplex.mat

Zi = projection(Zi,[1,2]);
ZiOrig = Zi;

G = Zi.G;
A = Zi.A;
b = Zi.b;
c = Zi.c;

tol = 1e-15;
tolNonZero = 1e-6;
tolR = 1e-15;
tolFeas = 1e-4;
tolE = 1e-15;

disp([num2str(Zi.nG) ' generators and ' num2str(Zi.nC) ' constraints (initial set).'])

% Zero columns and rows
zeroCols = find(all(abs([G;A])<tol,1));
G(:,zeroCols) = [];
A(:,zeroCols) = [];

zerorows = find(all(abs([A,b])<tol,2));
A(zerorows,:) = [];
b(zerorows) = [];

disp([num2str(Zi.nG) ' generators and ' num2str(Zi.nC) ' constraints (zeros removed).'])

% Remove factor as soon as redundancy detected
count = 0;
iterMax = 1000000; %100
for iter = 1:iterMax

    % Dimensions
    [nc, ng] = size(A); % nc = number of constraints, ng = number of generators

    generators_in=ng;

    % Initialize bounds
    E = repmat([-1, 1], ng, 1); % E initialized as [-1, 1] for all generators
    R = repmat([-inf, inf], ng, 1); % R initialized as [-Inf, Inf] for all generators

    % Iterate over constraints and generators
    for j = 1:ng       % Swapped the order of these for loops
        
        nonZeroIndx = find(abs(A(:,j))>=tolNonZero);

        % Vectorize the computation of R
        k = [1: j-1, j+1: ng];
        E_m = 1/2*(E(k,1)+E(k,2));
        E_r = 1/2*(E(k,2)-E(k,1));

        for i_indx = 1:length(nonZeroIndx)
            i = nonZeroIndx(i_indx);
            aij = A(i, j);

            % Vectorize the computation of R
            middleR = 1/aij*A(i,k)*E_m;
            rangeR = sum(abs(1/aij*A(i,k).*E_r'),2);

            % Calculate the bounds for constraint i
            Rj_lower = b(i) / aij - (middleR+rangeR);
            Rj_upper = b(i) / aij - (middleR-rangeR);

            % Update bounds for generator j
            R(j, 1) = max(R(j, 1), Rj_lower);
            R(j, 2) = min(R(j, 2), Rj_upper);

            redundant = (R(j,1) >= -1 + tolR) & (R(j,2) <= 1 - tolR);

            if redundant
                % Remove redundant generators and constraints simultaneously
                removable_constraint = i;
                remove_position = j;
                % Step 1: Define the E matrix
                E_lam = zeros(ng, nc);
                E_lam(remove_position, removable_constraint) = 1; % Place a 1 at the (remove_position, removable_constraint) position
                % Step 2: Compute \Lambda_G and \Lambda_A
                a_1j = A(removable_constraint, remove_position); % The removable constraint row
                Lambda_G = G * E_lam / a_1j; % \Lambda_G = G * E * a_1j^(-1)
                Lambda_A = A * E_lam / a_1j; % \Lambda_A = A * E * a_1j^(-1)
                % Step 3: Update the generator matrix G
                G = G - Lambda_G * A;
                % Step 4: Update the center vector c
                c = c + Lambda_G * b;
                % Step 5: Update the constraint matrix A
                A = A - Lambda_A * A;
                % Step 6: Update the constraint vector b
                b = b - Lambda_A * b;
                % Step 7: Remove the selected constraint from all matrices
                G(:, remove_position) = []; % Remove the column in G corresponding to remove_position
                A(removable_constraint, :) = []; % Remove the row in A corresponding to removable_constraint
                A(:, remove_position) = []; % Remove the column in A corresponding to remove_position
                b(removable_constraint,:) = [];
                break
            else
                % Update bounds for generator j in E
                E(j, 1) = max(E(j, 1), R(j, 1)); % Update lower bound
                E(j, 2) = min(E(j, 2), R(j, 2)); % Update upper bound

                % If E(j, 1) > E(j, 2), swap them
                if E(j, 1) > E(j, 2)
                    E(j, :) = E(j, [2, 1]); % Efficient swap in one step
                end
            end 
        end
        if redundant
            break
        end
    end

    % Zi = conZono(G,c,A,b);
    % disp([num2str(Zi.nG) ' generators and ' num2str(Zi.nC) ' constraints (redundancy removal).'])

    if redundant
        continue
    end

    E_m = 1/2*(E(:,1)+E(:,2));
    E_r = 1/2*(E(:,2)-E(:,1));
    E_r(E_r<=tolE) = 0;

    % Rescale
    c = c + G*E_m;
    G = G*diag(E_r);
    b = b - A*E_m;
    A = A*diag(E_r);


    % Zero columns and rows
    tol = 1e-10;
    zeroCols = find(all(abs([G;A])<tol,1));
    G(:,zeroCols) = [];
    A(:,zeroCols) = [];

    zerorows = find(all(abs([A,b])<tol,2));
    A(zerorows,:) = [];
    b(zerorows) = [];

    if isempty(zeroCols) && isempty(zerorows)

        R(:,1)=(R(:,1)-E_m)./E_r;
        R(:,2)=(R(:,2)-E_m)./E_r;
    
        totalExceed = max(0,-1-R(:,1))+max(0,-1+R(:,2));
        [totalExceedSorted,indxSorted] = sort(totalExceed);

        for m = 1:ng
            j=indxSorted(m);
            optSolver = 'mosek';  % Use Mosek solver
            s = zeros(ng, 1);  % Objective coefficients
            lbs = -ones(size(G,2), 1);  % Lower bounds for variables
            ubs = ones(size(G,2), 1);   % Upper bounds for variables
            
            % Set up the MOSEK optimization model
            prob.c = s;  % Objective coefficients (same as before)
            prob.a = sparse(A);  % Coefficient matrix
            prob.blc = b;  % Right-hand side of the constraints (lower bounds)
            prob.buc = b;  % Right-hand side of the constraints (upper bounds)
            
            % Feasibility checks (adjust bounds as needed)
            lbs(j) = 1 + tolFeas;  % Adjust the lower bound for the j-th variable
            ubs(j) = inf;  % Adjust the upper bound for the j-th variable
            prob.blx(j) = lbs(j);  % Update the lower bound for the j-th variable
            prob.bux(j) = ubs(j);  % Update the upper bound for the j-th variable
            
            % Set bounds on the variables
            prob.blx = lbs;  % Lower bounds on all variables
            prob.bux = ubs;  % Upper bounds on all variables
            
            % Perform the optimization
            [r, res] = mosekopt('maximize', prob);
            
            % Check the result status
            RjLPub = res.sol.bas.solsta;  % Get the solution status
            
            if strcmp(RjLPub, 'OPTIMAL')
                RjLPlb = RjLPub;  % If the solution is optimal, store the result status

                else
                   % Set up the MOSEK optimization model
                    prob.c = -s;  % Objective coefficients (same as before)
                    prob.a = sparse(A);  % Coefficient matrix
                    prob.blc = b;  % Right-hand side of the constraints (lower bounds)
                    prob.buc = b;  % Right-hand side of the constraints (upper bounds)
                    
                    % Feasibility checks (adjust bounds as needed)
                    lbs(j) = -inf;  % Adjust the lower bound for the j-th variable
                    ubs(j) = -1-tolFeas;  % Adjust the upper bound for the j-th variable
                    prob.blx(j) = lbs(j);  % Update the lower bound for the j-th variable
                    prob.bux(j) = ubs(j);  % Update the upper bound for the j-th variable
                    
                    % Set bounds on the variables
                    prob.blx = lbs;  % Lower bounds on all variables
                    prob.bux = ubs;  % Upper bounds on all variables
                    
                    % Perform the optimization
                    [r, res] = mosekopt('maximize', prob);
                    
                    % Check the result status
                    RjLPlb = res.sol.bas.solsta;  % Get the solution status
            end
                redundant = (~strcmp(RjLPub,'OPTIMAL')) && (~strcmp(RjLPlb,'OPTIMAL'));
            % end
  
            if redundant
                count = count + 1;
                % Remove redundant generators and constraints simultaneously
                removable_constraint = find(abs(A(:,j)) >= 1e-6,1);
                remove_position = j;
                % Step 1: Define the E matrix
                E_lam = zeros(ng, nc);
                E_lam(remove_position, removable_constraint) = 1; % Place a 1 at the (remove_position, removable_constraint) position
                % Step 2: Compute \Lambda_G and \Lambda_A
                a_1j = A(removable_constraint, remove_position); % The removable constraint row
                Lambda_G = G * E_lam / a_1j; % \Lambda_G = G * E * a_1j^(-1)
                Lambda_A = A * E_lam / a_1j; % \Lambda_A = A * E * a_1j^(-1)
                % Step 3: Update the generator matrix G
                G = G - Lambda_G * A;
                % Step 4: Update the center vector c
                c = c + Lambda_G * b;
                % Step 5: Update the constraint matrix A
                A = A - Lambda_A * A;
                % Step 6: Update the constraint vector b
                b = b - Lambda_A * b;
                % Step 7: Remove the selected constraint from all matrices
                G(:, remove_position) = []; % Remove the column in G corresponding to remove_position
                A(removable_constraint, :) = []; % Remove the row in A corresponding to removable_constraint
                A(:, remove_position) = []; % Remove the column in A corresponding to remove_position
                b(removable_constraint,:) = [];
                break
            end
        end
    end
generators_out=size(G,2);
if generators_in==generators_out
    break
end
end

Zi = conZono(G,c,A,b);
disp([num2str(Zi.nG) ' generators and ' num2str(Zi.nC) ' constraints (redundancy removal).'])

% Zi = conZono(G,c,A,b);
% figure;hold on
% tic;
% plot(ZiOrig,'r',0.1);
% toc
% tic;
% plot(Zi,'b',0.1)
% toc

