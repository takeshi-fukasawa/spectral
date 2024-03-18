function [x_sol_cell,other_output_k_2,iter_info]=...
    SQUAREM_func(fun_fp,spec,...
    x_0_cell,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% update_spec: dimension of a variable whose updating tune parameters  should be independently chosen (e.g. time)
%%% update_spec==0 => Implement standard fixed point iteration
% dampening_param

n_var=size(x_0_cell,2);
run preliminary_spectral.m

%%spec.compute_alpha_spec=1;

tic

% varargin:1*XXX


fun_k_cell={};

%% Read inputs
other_input_cell=varargin;


DIST_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);
    
%%%%%%%% Loop %%%%%%%%%%%
x_k_cell=x_0_cell;
count=0;

FLAG_ERROR=0;
for k=0:floor(ITER_MAX/2)-1

        [fun_k_1_cell,other_output_k_1]=...
           fun_fp(x_k_cell{:},other_input_cell{:});
        count=count+1;

       %%% Evaluate distance 
       for i=1:n_var
           r_k{1,i}=fun_k_1_cell{1,i}-x_k_cell{1,i};
           DIST_k_1(1,i)=norm_func(r_k{1,i}(:),spec.norm_spec);
           DIST_table(k*2+1,:)=DIST_k_1;

       end % for loop wrt i
       if sum(DIST_k_1(1,:))<TOL
            other_output_k_2=other_output_k_1;
            break;
       end % if statement

        [fun_k_2_cell,other_output_k_2]=...
           fun_fp(fun_k_1_cell{:},other_input_cell{:});
            count=count+1;

       for i=1:n_var
           v_k{1,i}=fun_k_2_cell{1,i}-2*fun_k_1_cell{1,i}+x_k_cell{1,i};

           difference=fun_k_2_cell{1,i}(:)-fun_k_1_cell{1,i}(:);
           DIST_k_2(1,i)=norm_func(difference(:),spec.norm_spec);
           
       end % for loop wrt i
       if sum(DIST_k_2(1,:))<TOL
           x_k_plus_1_cell=fun_k_2_cell;
           break;
       end % if statement

       if update_spec==0
          for i=1:n_var
              alpha_k{1,i}=1;% fixed point iterations
          end
      else
        alpha_k=compute_alpha_func(...
         r_k,v_k,...
            common_alpha_spec,compute_alpha_spec,dampening_param,update_spec);
     end

      for i=1:n_var
            x_k_plus_1_cell{1,i}=...
                x_k_cell{1,i}+2*alpha_k{1,i}.*r_k{1,i}+v_k{1,i}.*(alpha_k{1,i}).^2;
       end % for loop wrt i

    %%DIST_table(k+2,:)=DIST_k_2;
    DIST_table(k*2+2,:)=DIST_k_2;

    for i=1:n_var
        alpha_vec(1,i)=mean(alpha_k{i});
    end

    alpha_table(k+2,:)=alpha_vec;
    DIST=nanmax(DIST_k_2);

    if isnan(DIST)==1|isinf(sum(DIST))==1|isnan(sum(DIST))==1
       %warning("Error ?? ")
       x_k_plus_1_cell=x_k_cell;
       FLAG_ERROR=1;
       break;
    end
   
   interval=10;
    if k-floor(k/interval)*interval==0&DEBUG==1
        DIST_vec
    end

    if DIST<TOL
        %DIST
        FLAG_ERROR=0;
        break;
    end

    x_k_cell=x_k_plus_1_cell;

end %% end of for loop wrt k=0:ITER_MAX-1


%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc;
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=count;

iter_info.ITER_MAX=ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;

iter_info.DIST_table=DIST_table;
iter_info.alpha_table=alpha_table;


end


