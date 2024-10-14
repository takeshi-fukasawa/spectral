function vec=cell_to_vec_func(cell)

    n_vars=size(cell,2);
    vec=[];
    for i=1:n_vars
        vec=[vec;cell{i}(:)];
    end

end
