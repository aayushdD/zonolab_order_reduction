function out = pontry_sparse_diff_JK(obj1, obj2)
% Inputs:
%   X - struct with fields: c, G, A, b (constrained zonotope)
%   Y - struct with fields: G (zero-centered zonotope)
% Output:
%   result - projected constrained zonotope

out = obj1;
out.c = obj1.c - obj2.c;

for i = 1:obj2.nG
    out_plus = out;
    out_plus.c = out.c + obj2.G(:,i);
    out_minus = out;
    out_minus.c = out.c - obj2.G(:,i);
    out = sparse_intersection_JK(out_minus,out_plus);
end

% % Step 1: Construct Z1 and Z2 using the shifting centers method
% Z1.c = X.c - Y.G * ones(size(Y.G, 2), 1);
% Z1.G = X.G;
% Z1.A = X.A;
% Z1.b = X.b;
% 
% Z2.c = X.c + Y.G * ones(size(Y.G, 2), 1);
% Z2.G = X.G;
% Z2.A = X.A;
% Z2.b = X.b;
% 
% % Step 2: Compute the Cartesian product Z3 = Z1 × Z2
% Z3.c = [Z1.c; Z2.c];
% Z3.G = [Z1.G, zeros(size(Z1.G)); zeros(size(Z2.G)), Z2.G];
% Z3.A = blkdiag(Z1.A, Z2.A);
% Z3.b = [Z1.b; Z2.b];
% 
% % Step 3: Define Z4 as the zero constrained zonotope {0}
% Z4.c = zeros(2,1);
% Z4.G = zeros(2, 2); % 2 generators, all zero
% Z4.A = zeros(1, 2);
% Z4.b = 0;
% 
% % Step 4: Generalized intersection Z5 = Z3 ∩_{[I -I]} Z4
% % Construct the constraint matrices
% A1 = Z3.A;
% A2 = Z4.A;
% b1 = Z3.b;
% b2 = Z4.b;
% G1 = Z3.G;
% G2 = Z4.G;
% c1 = Z3.c;
% c2 = Z4.c;
% R = [eye(2), -eye(2)]; % [I -I]
% 
% Z5.c = c1;
% Z5.G = [G1, zeros(size(G1,1),size(G2,2))]; % G1 only (G2 = 0)
% Z5.A = [A1, zeros(size(A1,1),size(G2,2)); 
%         zeros(size(A2,1), size(G1,2)), A2;
%         R * G1, -G2];
% Z5.b = [b1; b2; c2 - R * c1];
% 
% % Step 5: Project Z5 using [I 0]
% dim = length(Z1.c);
% P = [eye(dim), zeros(dim)];
% Z_proj.c = P * Z5.c;
% Z_proj.G = P * Z5.G;
% Z_proj.A = Z5.A;
% Z_proj.b = Z5.b;
% 
% 
% result=conZono(Z_proj.G,Z_proj.c,Z_proj.A,Z_proj.b);
end

