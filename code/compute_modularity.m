function Q = compute_modularity(A, roi_list)

    N = size(A, 1);
  
    m = sum(A(:)) / 2;  
    

    degree = sum(A, 2);  

    Q = 0;
    

    for i = 1:N
        for j = 1:N
            if roi_list(i) == roi_list(j)  
                delta = 1;
            else
                delta = 0;
            end

            Q = Q + (A(i,j) - (degree(i) * degree(j)) / (2 * m)) * delta;
        end
    end
    

    Q = Q / (2 * m);
end
