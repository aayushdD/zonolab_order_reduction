function out = new_pontry_approach(X, Y)

    [V,~] = plot(Y);
    close all;

    for i = 1:size(V,1)
        v = V(i,:);
        Zv = zono(zeros(2,[]),-[v(1);v(2)]); 
        Xi = plus(X, Zv);
        if i == 1
            out = Xi;
        else
            out = sparse_intersection_JK(out , Xi);
        end
    end
end
 