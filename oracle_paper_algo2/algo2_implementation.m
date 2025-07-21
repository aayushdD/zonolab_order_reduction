 
 % load('pontry_example_reducable_example1.mat')
 % load('pontry_example_nonreducable_example1.mat')
 load('pontry_example_nonreducable_example2.mat')


 [v,~]=plot(Zi,'r',0.5);

 disp(Zi);

 [isZono, generators] = isZonotope_ds6(v);

 result_zono=zono(generators,[0;0]);

 disp(result_zono);

 hold on;
 plot(result_zono,'b',0.3);

 legend('actual','resulting');


