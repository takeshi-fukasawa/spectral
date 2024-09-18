function [x_sol_cell,other_output_k_2,iter_info]=...
    SQUAREM_func(fun_fp,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in May 2024.

n_var=size(x_0_cell,2);

spec=preliminary_spectral_func(spec,n_var);

ITER_MAX=spec.ITER_MAX;

%%spec.compute_alpha_spec=1;

t_SQUAREM=tic;

% varargin:1*XXX


fun_k_cell={};
FLAG_ERROR=0;

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
           DIST_k_1(1,i)=norm_func(r_k{1,i}(:),x_k_cell{1,i}(:),spec.norm_spec(i));
       end % for loop wrt i
       
       DIST_table(k*2+1,:)=DIST_k_1;

       if sum(DIST_k_1(1,:))<spec.TOL
            other_output_k_2=other_output_k_1;
            break;
       end % if statement

        [fun_k_2_cell,other_output_k_2]=...
           fun_fp(fun_k_1_cell{:},other_input_cell{:});
            count=count+1;

       for i=1:n_var
           v_k{1,i}=fun_k_2_cell{1,i}-2*fun_k_1_cell{1,i}+x_k_cell{1,i};

           difference=fun_k_2_cell{1,i}(:)-fun_k_1_cell{1,i}(:);
           DIST_k_2(1,i)=norm_func(difference(:),fun_k_1_cell{1,i}(:),spec.norm_spec(i));
           
       end % for loop wrt i
       if sum(DIST_k_2(1,:))<spec.TOL
           x_k_plus_1_cell=fun_k_2_cell;
           break;
       end % if statement

       if spec.update_spec==0
          for i=1:n_var
              alpha_k{1,i}=1;% fixed point iterations
          end
       else
        [alpha_k,alpha_max]=compute_alpha_func(...
         r_k,v_k,spec,k);
        end

      for i=1:n_var
            x_k_plus_1_cell{1,i}=...
                x_k_cell{1,i}+2*alpha_k{1,i}.*r_k{1,i}+v_k{1,i}.*(alpha_k{1,i}).^2;
       end % for loop wrt i

    %%DIST_table(k+2,:)=DIST_k_2;
    DIST_table(k*2+2,:)=DIST_k_2;

    for i=1:n_var
        %%alpha_vec(1,i)=mean(alpha_k{i}(:));
        alpha_vec(1,i)=max(alpha_k{i}(:));
        
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
    if k-floor(k/interval)*interval==0&spec.DEBUG==1
        DIST_k_2
    end

    if DIST<spec.TOL
        %DIST
        FLAG_ERROR=0;
        break;
    end

    x_k_cell=x_k_plus_1_cell;

end %% end of for loop wrt k=0:ITER_MAX-1


%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc(t_SQUAREM);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=count;

iter_info.ITER_MAX=ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;

iter_info.DIST_table=DIST_table;
iter_info.alpha_table=alpha_table;
iter_info.norm_spec=spec.norm_spec;


end


