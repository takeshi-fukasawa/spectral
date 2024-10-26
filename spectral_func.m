function [x_sol_cell,other_output_k_plus_1,iter_info]=...
    spectral_func(fun,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in May 2024.

n_var=size(x_0_cell,2);

spec=preliminary_spectral_func(spec,n_var);

other_input_cell=varargin;

FLAG_ERROR=0;

x_max_cell=spec.x_max_cell;
x_min_cell=spec.x_min_cell;

if spec.SQUAREM_spec==0

% varargin:1*XXX

%% Read inputs
ITER_MAX=spec.ITER_MAX;
DIST_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);
ITER_table_LINE_SEARCH=NaN(ITER_MAX,1);
step_size_table=NaN(ITER_MAX,n_var);

tic_spectral=tic;
 if spec.bound_spec==1
    x_0_cell=projection_func(x_0_cell,x_max_cell,x_min_cell);
 end
 [fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});

 if spec.fixed_point_iter_spec==1
     for i=1:n_var
         fun_0_cell{1,i}=fun_0_cell{1,i}-x_0_cell{1,i};
     end
 end


feval=1;
    
    %%% DIST: sup norm of F(x)=x-Phi(x). 
    DIST_vec=ones(1,n_var);
    for i=1:n_var
      DIST_vec(1,i)=norm_func(fun_0_cell{1,i}(:),x_0_cell{1,i}(:),spec.norm_spec(i));

      alpha_0{1,i}=spec.alpha_0;
      
      if spec.update_spec==0
        alpha_0{1,i}=1;
      end

      alpha_table(1,i)=alpha_0{1,i};
      
    end % loop wrt i


    DIST=nanmax(DIST_vec);
    DIST_table(1,:)=DIST_vec;


    obj_val_table=NaN(ITER_MAX,n_var);
    obj_val_table(1,:)=DIST_vec.^2;% L2 norm

conv=(sum((DIST_vec<spec.TOL),'all')==n_var);

x_k_cell=x_0_cell;
fun_k_cell=fun_0_cell;
other_output_k_plus_1=other_output_0;

x_k_minus_1_cell=x_0_cell;%%%
fun_k_minus_1_cell=fun_0_cell;%%%

%%%%%%%% Loop %%%%%%%%%%%

if conv==0 & ITER_MAX>=2
for k=0:ITER_MAX-1


   if k>=1
      for i=1:n_var
        Delta_x_cell{1,i}=x_k_cell{i}-x_k_minus_1_cell{i};
        Delta_fun_cell{1,i}=fun_k_cell{i}-fun_k_minus_1_cell{i};
      end % loop wrt i


      if spec.update_spec==0
          for i=1:n_var
              alpha_k{1,i}=1;
          end
      else
        [alpha_k,alpha_max]=compute_alpha_func(...
         Delta_x_cell,Delta_fun_cell,spec,k);
         spec.alpha_max=alpha_max;
      end

    %%% Update variables %%%%%%%%%%%%%%%

    for i=1:n_var
        d_k_cell{1,i}=alpha_k{1,i}.*fun_k_cell{1,i};
    end % for loop wrt i



    [x_k_plus_1_cell, fun_k_plus_1_cell,...
    other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec,...
    obj_val_vec,step_size]=...
        spectral_update_func(fun,x_k_cell,fun_k_cell,alpha_k,d_k_cell,other_input_cell,...
        n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);

    ITER_table_LINE_SEARCH(k+2,1)=iter_line_search;%% Number of line search iterations

    feval=feval+iter_line_search;

    DIST_table(k+2,:)=DIST_vec;
    alpha_table(k+2,:)=alpha_vec;
    obj_val_table(k+2,:)=obj_val_vec;
    step_size_table(k+2,:)=step_size;
    DIST=nanmax(DIST_vec);
    

    if isnan(DIST)==1|isinf(sum(DIST_table(k+2,:)))==1|isnan(sum(DIST_table(k+2,:)))==1
       %warning("Error ?? ")
       x_k_plus_1_cell=x_k_cell;
       FLAG_ERROR=1;
       break;
    end
   
   interval=10;
    if k-floor(k/interval)*interval==0&spec.DEBUG==1
        DIST_vec
    end

    if sum((DIST_vec<spec.TOL),'all')==n_var
        FLAG_ERROR=0;
        %DIST
        break;
    end

    %%% Replace variables for the next iteration
    x_k_minus_1_cell=x_k_cell;
    fun_k_minus_1_cell=fun_k_cell;


	x_k_cell=x_k_plus_1_cell;
    fun_k_cell=fun_k_plus_1_cell;
   

    
  else % k==0
    for i=1:n_var
     alpha_k{1,i}=alpha_0{1,i};
     if isempty(spec.dampening_param)==0
      alpha_k{1,i}=alpha_k{1,i}*(spec.dampening_param{1,i});
     end


   end% for loop wrt i

   
 end%k>=1 or k==0
end %% end of for loop wrt k=0:ITER_MAX-1


else % no iteration
    k=0;
    x_k_plus_1_cell=x_0_cell;
    other_output_k_plus_1=other_output_0;
    fun_k_cell=fun_0_cell;
end

%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc(tic_spectral);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=feval;
iter_info.ITER_MAX=spec.ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;
iter_info.obj_val_table=obj_val_table;

iter_info.fun_cell=fun_k_cell;

iter_info.DIST_table=DIST_table;
iter_info.alpha_table=alpha_table;
iter_info.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;
iter_info.norm_spec=spec.norm_spec;

iter_info.step_size_table=step_size_table;
iter_info.spec=spec;

elseif spec.SQUAREM_spec==1
    [x_sol_cell,other_output_k_plus_1,iter_info]=...
    SQUAREM_func(fun,spec,x_0_cell,other_input_cell{:});

elseif spec.Anderson_acceleration==1 % Anderson
    [x_sol_cell,other_output_k_plus_1,iter_info]=...
    Anderson_func(fun,spec,x_0_cell,other_input_cell{:});
    % Anderson code based on Zhang et al. (2020)
    % [x_sol_cell,other_output_k_plus_1,iter_info] = Anderson_acceleration_func(x_0_cell, fun,spec,...
    %     other_input_cell{:});
     
end



end

