%running script for new pontry approach

% X=zono(eye(2)*10,[0;0]);
% Y=zono(eye(2)*2,[0;0]);

load zono1.mat
load zono2.mat

zonolab_pontry=pontryDiff(Z1,Z2);
new_pontry=new_pontry_approach(Z1,Z2);
close all

figure;
hold on;
plot(Z1,'g',0.5);
plot(Z2,'k');
plot(zonolab_pontry,'r',0.5);
plot(new_pontry,'b',0.5);

legend('Z1','Z2','zonolab','new');

figure;
imagesc(new_pontry.A);
figure;
imagesc(zonolab_pontry.A);
%%
% G = [1 0; 0 1];
G = [1, 0.5, -0.5; 
      0, 0.9, 0.9];

c=zeros(2,1);
Z1=zono(G,c);


G = [0.2  0.1 -0.1  0.15 -0.15;  
      0.1 -0.2  0.15 0.1  -0.1];

% G = [0.1, 0.15, -0.05; 
%       0.05, -0.1, 0.2];

c=zeros(2,1);
Z2=zono(G,c);

Zi = pontryDiff(Z1,Z2);
Zj=new_pontry_approach(Z1,Z2);


figure;
hold on;
plot(Z1,'g',0.5);
plot(Z2,'k');
plot(Zi,'r',0.8);
plot(Zj,'b',0.5);

legend('Z1','Z2','zonolab','new');

figure;
imagesc(Zj.A);
figure;
imagesc(Zi.A);