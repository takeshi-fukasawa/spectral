function x_cell=...
    projection_func(x_cell,x_max_cell,x_min_cell)
 %% Restrict the range of x
   n_var=size(x_cell,2);
 
   for i=1:n_var
        if isempty(x_max_cell{1,i})==0
             x_cell{1,i}=min(x_cell{1,i},x_max_cell{1,i});

        end %isempty(x_max_cell{1,i})==0

        if isempty(x_min_cell{1,i})==0
             x_cell{1,i}=max(x_cell{1,i},x_min_cell{1,i});
       end%isempty(x_min_cell{1,i})==0
   end% for loop

end

