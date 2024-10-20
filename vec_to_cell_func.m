function x_cell=vec_to_cell_func(x_vec,elem_x,x_0_cell)

n_var=size(elem_x(:),1);
    loc=1;
    for i=1:n_var
        x_cell{1,i}=reshape(x_vec(loc:loc+elem_x(1,i)-1),size(x_0_cell{i}));
        loc=loc+elem_x(1,i);
    end

end
