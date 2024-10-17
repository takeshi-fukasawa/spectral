function spec=preliminary_spectral_func(spec,n_var);

if isfield(spec,'x_max_cell')==0
    spec.x_max_cell=cell(1,n_var);
end

if isfield(spec,'x_min_cell')==0
    spec.x_min_cell=cell(1,n_var);
end

if isfield(spec,'fixed_point_iter_spec')==0
    spec.fixed_point_iter_spec=1;%%%%%
end

if isfield(spec,'update_spec')==0
    spec.update_spec=[];
end

if isfield(spec,'dim_hetero_alpha')==0
    spec.dim_hetero_alpha=[];
elseif size(spec.dim_hetero_alpha(:),1)==1
    spec.dim_hetero_alpha=ones(n_var,1)*(spec.dim_hetero_alpha);
end

if isfield(spec,'dampening_param')==0
    spec.dampening_param=[];
end


if isfield(spec,'alpha_0')==0
    spec.alpha_0=1;
end


if isfield(spec,'common_alpha_spec')==0
    spec.common_alpha_spec=0;
end

if isfield(spec,'TOL')==0
    spec.TOL=1e-10;
end


if isfield(spec,'ITER_MAX')==0
    spec.ITER_MAX=1000;
end

if isfield(spec,'line_search_spec')==0
    spec.line_search_spec=0;
end

if isfield(spec,'ITER_MAX_LINE_SEARCH')==0
    spec.ITER_MAX_LINE_SEARCH=10;
end

if spec.line_search_spec==0
    spec.ITER_MAX_LINE_SEARCH=1;
end


if isfield(spec,'DEBUG')==0
    spec.DEBUG=0;
end

if isfield(spec,'compute_alpha_spec')==0
    spec.compute_alpha_spec=3;
end

if isfield(spec,'SQUAREM_spec')==0
    spec.SQUAREM_spec=0;%%%
end

if isfield(spec,'norm_spec')==0
    spec.norm_spec=0;
end

if spec.line_search_spec==1
    spec.norm_spec=2;
end

if size(spec.norm_spec(:),1)==1
    spec.norm_spec=ones(1,n_var).*(spec.norm_spec);
end


spec.bound_spec=0;
for i=1:n_var
   if isempty(spec.x_max_cell{1,i})==0 | isempty(spec.x_min_cell{1,i})==0
     spec.bound_spec=1;
   end
end

if isfield(spec,'positive_alpha_spec')==0
    spec.positive_alpha_spec=0;%%%%%%
end


if isfield(spec,'alpha_max')==0
    spec.alpha_max=10^10;
end

if isfield(spec,'alpha_min')==0
    spec.alpha_min=-10^10;
end


if spec.positive_alpha_spec==1
    if isfield(spec,'alpha_min')==0
        spec.alpha_min=1e-8; %% Restrict alpha to be positive
    end
    spec.compute_alpha_spec=3;
end

%%%% Line search parameters
if isfield(spec,'M')==0
    spec.M=10;
end

if isfield(spec,'gamma')==0
    spec.gamma=10^(-4);
end

if isfield(spec,'rho')==0
    spec.rho=0.3;
end

if isfield(spec,'Anderson_acceleration')==0
    spec.Anderson_acceleration=0;
end

if spec.Anderson_acceleration==1
    spec.SQUAREM_spec=[];
end

end % function



