

% G = [1 0; 0 1];
G = [1, 0.5, -0.5; 
      0, 0.9, 0.9];

c=zeros(2,1);
Z1=zono(G,c);


% G = [0.2  0.1 -0.1  0.15 -0.15;  
%       0.1 -0.2  0.15 0.1  -0.1];

G = [0.1, 0.15, -0.05; 
      0.05, -0.1, 0.2];

c=zeros(2,1);
Z2=zono(G,c);

Zi = new_pontry_approach(Z1,Z2);