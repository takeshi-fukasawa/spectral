function [x_sol_cell,other_output_k_plus_1,conv_table,iter_info,fun_k_cell]=...
    spectral_bound_func(fun,spec,...
    x_0_cell,x_max_cell,x_min_cell,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% update_spec: dimension of a variable whose updating tune parameters  should be independently chosen (e.g. time)
%%% update_spec==0 => Implement standard fixed point iteration
% dampening_param

n_var=size(x_0_cell,2);

run preliminary_spectral.m


% varargin:1*XXX

fun_k_cell={};

%% Read inputs
other_input_cell=varargin;

bound_spec=0;
for i=1:n_var
   if isempty(x_max_cell{1,i})==0 | isempty(x_min_cell{1,i})==0
     bound_spec=1;
   end
end


DIST_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);
ITER_table_LINE_SEARCH=NaN(ITER_MAX,1);

FLAG_ERROR=[];


tic
 if bound_spec==0         
     [fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
else
    [x_0_cell,fun_0_cell,other_output_0]=fun_bdd_func(...
        fun,x_0_cell,other_input_cell,x_max_cell,x_min_cell,n_var,max_opt_spec);
end

    
    %%% DIST: sup norm of F(x)=x-Phi(x). 
    DIST_vec=ones(1,n_var);
    for i=1:n_var
      DIST_vec(1,i)=max(abs(fun_0_cell{1,i}),[],'all','omitnan');
    end

    DIST=nanmax(DIST_vec);
    DIST_table(1,:)=DIST_vec;
    alpha_table(1,:)=alpha_0;

    

x_k_cell=x_0_cell;
fun_k_cell=fun_0_cell;
other_output_k_plus_1=other_output_0;

%%%%%%%% Loop %%%%%%%%%%%

if DIST>TOL
for k=0:ITER_MAX-1

   %%% In the current iteration, fun_k_cell is assumed to be far from zero

    for i=1:n_var
     if k>=1
     Delta_x_i=x_k_cell{i}-x_k_minus_1_cell{i};
     Delta_fun_i=fun_k_cell{i}-fun_k_minus_1_cell{i};

      if isempty(update_spec)==1
        sum_dim_ids='all';

      elseif sum(update_spec(:))>0 %%% XXX-dependent tune parameters
         sum_dim_ids=1:size(size(Delta_x_i),2);
         sum_dim_ids=sum_dim_ids(sum_dim_ids~=update_spec(i));
      end

      if isempty(update_spec)==1 | (isempty(update_spec)==0 & sum(update_spec(:))>0) 
        sum_Delta_x_x_i=sum(Delta_x_i.*Delta_x_i,sum_dim_ids,'omitnan');%vector
        sum_Delta_fun_fun_i=sum(Delta_fun_i.*Delta_fun_i,sum_dim_ids,'omitnan');%vector      
        sum_Delta_x_fun_i=sum(Delta_x_i.*Delta_fun_i,sum_dim_ids,'omitnan');%vector  
      end

     if common_alpha_spec==1
      sum_Delta_x_fun_cell{1,i}=sum_Delta_x_fun_i;
      sum_Delta_x_x_cell{1,i}=sum_Delta_x_x_i;
      sum_Delta_fun_fun_cell{1,i}=sum_Delta_fun_fun_i;
     end

     if isempty(update_spec)==1 | (isempty(update_spec)==0 & sum(update_spec(:))>0) %%%
        
        alpha_k_i=compute_alpha_func(sum_Delta_x_fun_i,...
            sum_Delta_x_x_i,sum_Delta_fun_fun_i,...
            compute_alpha_spec,dampening_param);

     else
        alpha_k_i=1;
     end

   else%% k==1
     alpha_k_i=alpha_0; 
   end
    
      alpha_k{1,i}=alpha_k_i;

   end % for loop wrt i

   if common_alpha_spec==1 & k>=2
       sum_Delta_x_fun=sum([sum_Delta_x_fun_cell{:}]);
       sum_Delta_x_x=sum([sum_Delta_x_x_cell{:}]);
       sum_Delta_fun_fun=sum([sum_Delta_fun_fun_cell{:}]);

       alpha_k_val=compute_alpha_func(sum_Delta_x_fun,...
       sum_Delta_x_x,sum_Delta_fun_fun,...
       compute_alpha_spec,dampening_param);

       for i=1:n_var
           alpha_k{1,i}=alpha_k_val;
       end

   end %common_alpha_spec==1 & k>=2
   

     
    %%% Update variables %%%%%%%%%%%%%%%
    [x_k_plus_1_cell, fun_k_plus_1_cell,...
    other_output_k_plus_1,DIST_vec,iter_line_search,alpha_vec]=...
        spectral_update_func(fun,x_k_cell,alpha_k,fun_k_cell,other_input_cell,...
        n_var,line_search_spec,...
        DIST_table,ITER_MAX_LINE_SEARCH,bound_spec,...
        x_max_cell,x_min_cell,k,max_opt_spec);

    ITER_table_LINE_SEARCH(k+2,1)=iter_line_search;%% Number of line search iterations

    DIST_table(k+2,:)=DIST_vec;
    alpha_table(k+2,:)=alpha_vec;
    DIST=nanmax(DIST_vec);

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

    %%% Replace variables for the next iteration
	x_k_minus_1_cell=x_k_cell;
	x_k_cell=x_k_plus_1_cell;
	fun_k_minus_1_cell=fun_k_cell;
    fun_k_cell=fun_k_plus_1_cell;
   
end %% end of for loop wrt k=0:ITER_MAX-1


else % no iteration
    k=0;
    x_k_plus_1_cell=x_k_cell;
    x_sol_cell=x_0_cell;
    other_output_k_plus_1=other_output_0;
    DIST_MAT=[];
    fun_k_cell=fun_0_cell;
end

%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc;
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k;
iter_info.ITER_MAX=ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;

conv_table.DIST_table=DIST_table;
conv_table.alpha_table=alpha_table;
conv_table.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;

return


