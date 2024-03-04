function [x_sol_cell,other_output_k_plus_1,conv_table,iter_info,fun_k_cell]=...
    spectral_bound_func(fun,spec,...
    x_0_cell,x_max_cell,x_min_cell,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% vec: dimension of a variable whose updating tune parameters  should be independently chosen (e.g. time)
%%% vec==0 => Implement standard fixed point iteration
% dampening_param

tic

global FLAG_ERROR DIST count


%%%spec=[];%%%%%
run preliminary_spectral.m

n_var=size(x_0_cell,2);

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


 if bound_spec==0         
     [fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
else
    [x_0_cell,fun_0_cell,other_output_0]=fun_bdd_func(...
        fun,x_0_cell,other_input_cell,x_max_cell,x_min_cell,n_var);
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

eps_val=1e-6;
eps_val=0;


if DIST>TOL
for k=0:ITER_MAX-1

   %%% In the current iteration, fun_k_cell is assumed to be far from zero

    for i=1:n_var
     if k>=1
     Delta_x_i=x_k_cell{i}-x_k_minus_1_cell{i};
     Delta_fun_i=fun_k_cell{i}-fun_k_minus_1_cell{i};

      if isempty(vec)==1
        sum_dim_ids='all';

      elseif sum(vec(:))>0 %%% XXX-dependent tune parameters
         sum_dim_ids=1:size(size(Delta_x_i),2);
         sum_dim_ids=sum_dim_ids(sum_dim_ids~=vec(i));
      end

      if isempty(vec)==1 | (isempty(vec)==0 & sum(vec(:))>0) 
        sum_Delta_x_x=sum(Delta_x_i.*Delta_x_i,sum_dim_ids,'omitnan');%vector
        sum_Delta_fun_fun=sum(Delta_fun_i.*Delta_fun_i,sum_dim_ids,'omitnan');%vector      
        sum_Delta_x_fun=sum(Delta_x_i.*Delta_fun_i,sum_dim_ids,'omitnan');%vector
      
        %%sign_i=sign(sum_Delta_x_fun);%1*1*n_dim etc.
        numer_i=sqrt(sum_Delta_x_x)+eps_val;
        denom_i=sqrt(sum_Delta_fun_fun)+eps_val;
      end

     if common_alpha_spec==1
      sum_Delta_x_fun_cell{1,i}=sum_Delta_x_fun;
      sum_Delta_x_x_cell{1,i}=sum_Delta_x_x;
      sum_Delta_fun_fun_cell{1,i}=sum_Delta_fun_fun;
     end

     if isempty(vec)==1 | (isempty(vec)==0 & sum(vec(:))>0) %%%
        %%sign_i=1;%%%%%% same sign spec %%%
        %%alpha_k_i=sign_i.*numer_i./denom_i;%scalar or vector (wrt the dimension specified in "vec")
        alpha_k_i=numer_i./denom_i;%scalar or vector (wrt the dimension specified in "vec")
     else
        alpha_k_i=1;
     end

    alpha_k_i((isnan(alpha_k_i)==1))=1;%%%
    alpha_k_i((isinf(alpha_k_i)==1))=1;%%%
    alpha_k_i(((alpha_k_i==0)))=1;%%%

   else%% k==1
     alpha_k_i=alpha_0; 
   end

    %%alpha_k_i=0.1;
    
    if isempty(dampening_param)==0
        alpha_k_i=alpha_k_i*dampening_param(i);
    end

      alpha_k{1,i}=alpha_k_i;

   end % for loop wrt i

   if common_alpha_spec==1 & k>=2
       sum_Delta_x_fun=sum([sum_Delta_x_fun_cell{:}]);
       sum_Delta_x_x=sum([sum_Delta_x_x_cell{:}]);
       sum_Delta_fun_fun=sum([sum_Delta_fun_fun_cell{:}]);

       numer=sqrt(sum_Delta_x_x)+eps_val;
       denom=sqrt(sum_Delta_fun_fun)+eps_val;
       alpha_k_val=numer./denom;

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
        x_max_cell,x_min_cell,k);

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
count=k;
%DIST_vec

else % no iteration
    warning("XXX")
    count=1;
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

conv_table.DIST_table=DIST_table;
conv_table.alpha_table=alpha_table;
conv_table.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;

return


