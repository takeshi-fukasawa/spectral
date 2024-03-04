function [x_sol_cell,other_output_k_plus_1,conv_table,iter_info,fun_k_cell]=...
    SQUAREM_func(fun,n_var,vec,dampening_param,...
    x_0_cell,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% vec: dimension of a variable whose updating tune parameters  should be independently chosen (e.g. time)
%%% vec==0 => Implement standard fixed point iteration
% dampening_param
%%% function: fixed point iteration representation

tic

global DEBUG FLAG_ERROR DIST count ITER_MAX TOL

TOL=1e-12;
common_alpha_spec=0;

ITER_MAX=3000;

% varargin:1*XXX

fun_k_cell={};

%% Read inputs
other_input_cell=varargin;


DIST_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);
    
%%%%%%%% Loop %%%%%%%%%%%
x_k_cell=x_0_cell;

eps_val=1e-6;
eps_val=0;


if DIST>TOL
for k=0:ITER_MAX-1

        [fun_k_1_cell,other_output_k_1]=...
           fun(x_k_cell{:},other_input_cell{:});

       %%% Evaluate distance 
       for i=1:n_var
           r_k{1,i}=fun_k_1_cell{1,i}-x_k_cell{1,i};
           DIST_k_1(1,i)=max(abs(r_k{1,i}(:)));
       end % for loop wrt i
       if sum([DIST_k_1_cell{1,i}{:}])<TOL
        break;
       end

        [fun_k_2_cell,other_output_k_2]=...
           fun(x_k_1_cell{:},other_input_cell{:});

       for i=1:n_var
           v_k{1,i}=fun_k_2_cell{1,i}-2*fun_k_2_cell{1,i}+x_k_cell{1,2};
           DIST_k_2(1,i)=max(abs(v_k{1,i}(:)));
       end % for loop wrt i
       if sum([DIST_k_1_cell{1,i}{:}])<TOL
        break;
       end

   %%% alpha=-1 ??? %%%%%

   %%% Compute alpha
    for i*1:n_var
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

end %% end of for loop wrt k=0:ITER_MAX-1
count=k;
%DIST_vec

end

%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc;
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k;
iter_info.ITER_MAX=ITER_MAX;

conv_table.DIST_table=DIST_table;
conv_table.alpha_table=alpha_table;

return


