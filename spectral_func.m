function [x_sol_cell,other_output_k_plus_1,iter_info]=...
    spectral_func(fun,spec,...
    x_0_cell,varargin)

%%% Written by Takeshi Fukasawa in May 2024.

n_var=size(x_0_cell,2);

run preliminary_spectral.m

other_input_cell=varargin;

FLAG_ERROR=0;

if spec.SQUAREM_spec==0

% varargin:1*XXX

%% Read inputs

DIST_table=NaN(ITER_MAX,n_var);
obj_val_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);
ITER_table_LINE_SEARCH=NaN(ITER_MAX,1);
step_size_table=NaN(ITER_MAX,n_var);

tic_spectral=tic;
if spec.bound_spec==1
    x_0_cell=projection_func(x_0_cell,x_max_cell,x_min_cell);
end

if spec.with_obj_val_spec==0
    [fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
else
    [obj_fun_0_cell,fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
end

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

    if isfield(spec,'alpha_0_spec')==0
        alpha_0{1,i}=1;
    elseif isfield(spec,'alpha_0_spec')==1
        if spec.alpha_0_spec==2
            max_val_i(1,i)=max(abs(fun_0_cell{1,i}(:)));
        else
            alpha_0{1,i}=1;
        end
    else
        alpha_0{1,i}=1;
    end


    if spec.update_spec==0
        alpha_0{1,i}=1;
    end

    if isempty(alpha_0_param)==0
        alpha_0{1,i}=alpha_0_param;
    end
    alpha_table(1,i)=alpha_0{1,i};
      
end % loop wrt i

    DIST_table(1,:)=DIST_vec;
    if spec.with_obj_val_spec==0
        obj_val_table(1,:)=DIST_vec.^2;% L2 norm
    else%with_obj_val_spec==1
        for i=1:n_var
            obj_val_table(1,i)=obj_fun_0_cell{i};
        end
    end

    conv=(sum((DIST_vec<TOL),'all')==n_var);

x_k_cell=x_0_cell;
fun_k_cell=fun_0_cell;
other_output_k_plus_1=other_output_0;



%%%%%%%% Loop %%%%%%%%%%%

if conv==0 & ITER_MAX>=2
for k=0:ITER_MAX-2


   if k>=1

      for i=1:n_var
        Delta_x_cell{1,i}=x_k_cell{i}-x_k_minus_1_cell{i};
        Delta_fun_cell{1,i}=fun_k_cell{i}-fun_k_minus_1_cell{i};
      end % loop wrt i

      if spec.update_spec==0
          for i=1:n_var
              alpha_k{1,i}=1;
          end

      elseif spec.BFGS_spec==0
        [alpha_k,alpha_max]=compute_alpha_func(...
         Delta_x_cell,Delta_fun_cell,spec,k,DIST_table(k+1,:));
         spec.alpha_max=alpha_max;

      elseif spec.BFGS_spec==1
         s=Delta_x_cell{1,1}(:);

         y=Delta_fun_cell{1,1}(:);
         
         if spec.fixed_point_iter_spec==1
            y=y.*(-1);%%%####
         end

         I=eye(size(s,1));
         syp=s*y';
         spy=s'*y;% scalar
         ysp=y*s';
         ssp=s*s';

         H_k=(I-syp./spy)*H_k_minus_1*(I-ysp./spy)+ssp./spy;

         alpha_k{1,1}=1;%%%%%

    
    end

  else % k==0
      for i=1:n_var
       alpha_k{1,i}=alpha_0{1,i};
     end% for loop wrt i


      if spec.BFGS_spec==1
         H_k=eye(size(fun_k_cell{1,i}(:),1));
      end
   end

    if isempty(spec.dampening_param)==0
       dampening_param=spec.dampening_param;
       for i=1:n_var
        alpha_k{1,i}=alpha_k{1,i}*dampening_param{1,i};
       end
    end


    %%% Update variables %%%%%%%%%%%%%%%
   for i=1:n_var
        %%alpha_k_i{1,i}=1;%%%%%%%%%%%

        d_k_cell{1,i}=alpha_k{1,i}.*fun_k_cell{1,i};

        if spec.minimization_spec(i)==1
            d_k_cell{1,i}=-d_k_cell{1,i};% Direction of descending
        end

        
        if spec.conjugate_gradient_spec==1 & k>=1
            d_k_cell{1,i}=conjugate_gradient_func(...
            d_k_cell,i,alpha_k,Delta_fun_cell,Delta_x_cell,fun_k_cell,fun_k_minus_1_cell);
        end

        if spec.BFGS_spec==1 & k>=2
            %H_k=eye(size(d_k_cell{1,i}(:),1));
            
            d_k_cell{1,i}=H_k*d_k_cell{1,i}(:);
        end

        alpha_table(k+1,i)=max(alpha_k{1,i}(:));
    end % for loop wrt i

    [x_k_plus_1_cell, fun_k_plus_1_cell,...
    other_output_k_plus_1,DIST_vec,obj_val_vec,iter_line_search,step_size]=...
        spectral_update_func(fun,x_k_cell,fun_k_cell,d_k_cell,other_input_cell,...
        n_var,spec,x_max_cell,x_min_cell,k,obj_val_table);

        
    ITER_table_LINE_SEARCH(k+2,1)=iter_line_search;%% Number of line search iterations

    feval=feval+iter_line_search;



    DIST_table(k+2,:)=DIST_vec;
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
    if k-floor(k/interval)*interval==0&DEBUG==1
        DIST_vec

        if spec.with_obj_val_spec==1
            obj_val_vec
        end
    end

    if sum((DIST_vec<TOL),'all')==n_var
        FLAG_ERROR=0;
        %DIST
        break;
    end

    %%% Replace variables for the next iteration
	x_k_minus_1_cell=x_k_cell;
	x_k_cell=x_k_plus_1_cell;
	fun_k_minus_1_cell=fun_k_cell;
    fun_k_cell=fun_k_plus_1_cell;
   

    if spec.BFGS_spec==1
       H_k_minus_1=H_k;
    end


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

t_cpu=toc(tic_spectral);
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k+1;
iter_info.feval=feval;
iter_info.ITER_MAX=ITER_MAX;
iter_info.FLAG_ERROR=FLAG_ERROR;

iter_info.fun_cell=fun_k_cell;

iter_info.DIST_table=DIST_table;
iter_info.alpha_table=alpha_table;
iter_info.ITER_table_LINE_SEARCH=ITER_table_LINE_SEARCH;
iter_info.step_size_table=step_size_table;
iter_info.norm_spec=spec.norm_spec;

if spec.with_obj_val_spec==1
    iter_info.obj_val_table=obj_val_table;
end

else % spec.SQUAREM_spec==1

    [x_sol_cell,other_output_k_plus_1,iter_info]=...
    SQUAREM_func(fun,spec,x_0_cell,other_input_cell{:});
end

return


