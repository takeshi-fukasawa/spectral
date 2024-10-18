# MATLAB code of the Spectral, SQUAREM, and the Anderson acceleration algorithm
This repository contains the MATLAB code of the spectral algorithm (e.g., Barzilai and Borwein, 1988; La Cruz et al., 2006), the SQUAREM algorithm (e.g., Varadhan and Roland, 2008), and the Anderson acceleration method (e.g., Anderson, 1995) for solving nonlinear equations and fixed point problems, based on the discussion in Section 5 of Fukasawa (2024).

Spectral_func is the main function. It is designed to apply the spectral algorithm to solve nonlinear equations or fixed point problems.

It also introduces variable-type-specific step sizes in the spectral algorithm, which sometimes accelerate the convergence. For details, see Section 5 of Fukasawa (2024).

Besides, if we let `spec.SQUAREM_spec==1`, we can alternatively use the SQUAREM algorithm (SQUAREM_func). If we let `spec.Anderson_acceleration==1`, we can apply the Anderson acceleration (Anderson_func). 

The section "Examples" concisely describes the usage of the function. My other repository takeshi-fukasawa/BLP_algorithm shows the usage of the function, in the context of BLP estimation (estimation of random coefficient logit models using aggregate data) in economics. 

### Input of spectral_func
The essential input of the function is as follows:
* fun: fixed-point mapping we want to solve (fun(x_sol)=x_sol), or nonlinear function that should satisfy fun(x_sol)=0.
Note that there should be two outputs regarding the function. The first is the cell of variables, and the second one is a structure array containing other variables.
The second is the middle output of the function, and it can be used for further analysis. If there is no other variables, let the second output be []. For details of the usage, see also examples shown below.
* spec (Structure array): Detailed specifications of the spectral algorithm. See "Detailed specification of the spectral algorithm" for details.
If there is no problem in choosing the default values, setting spec=[] is sufficient.
* x_0_cell (Cell):  Initial values of variables

### Output of spectral_func
The output of the function is:
* output_cell (cell): Variables that should be solved
* other_vars: Other variables, which are not in output_cell
* iter_info (structure array): Information on the iteration;

  + t_cpu: Computation time
  + n_iter: Number of iterations
  + feval: Number of function evaluations
  + ITER_MAX: Maximum number of iterations specified in the function
  + FLAG_ERROR: If the algorithm converges or diverges, FLAG_ERROR=1. Otherwise, FLAG_ERROR=0.
  + fun_cell(cell): the values of fun(x)
  + DIST_table (ITER_MAX by n_var matrix):  the values of DIST (used to determine the convergence of the algorithm) in the iteration



## Examples
The following examples are from my other repository takeshi-fukasawa/BLP_algorithm, in the context of BLP estimations in economics.

### Example 1 (The case with only one type of variable)
Suppose we want to solve a fixed point problem concerning delta, represented by the following code in MATLAB:

```
[delta,other_vars]=BLP_update_func(delta,...
    weight,mu_ijt_est,rho_est,S_jt_data,tune_param_BLP);
```

Here, (weight,mu_ijt_est,rho_est,S_jt_data,tune_param_BLP) are exogenous additional parameters.

Then, we can numerically solve for delta by the following code:
```
[output_BLP_spectral,other_output,iter_info]=...
    spectral_func(@BLP_update_func,spec,{delta_initial0},...
    weight,mu_ijt_est,rho_est,S_jt_data,tune_param_BLP);

delta_sol=output_BLP_spectral{1};
```

Here, delta_initial0 denotes the initial values of delta, which is an array.

### Example 2 (The case with several types of variables)
Suppose we want to solve a fixed point problem concerning two types of variables delta and V, represented by the following code in MATLAB:

```
[out,other_vars]=BLP_Bellman_joint_update_func(delta,V,...
weight,mu_ijt_est,rho_est,...
    S_jt_data,weight_V,x_V,beta_C,L,tune_param_BLP);
delta=out{1};
V=out{2};
```

Then, we can numerically  solve for delta and V by the following code:
```
[output_spectral,other_vars,iter_info]=...
        spectral_func(@BLP_Bellman_joint_update_func,spec,...
        {delta_initial0,V_initial0},...
        weight,mu_ijt_est,rho_est,...
    S_jt_data,weight_V,x_V,beta_C,L,tune_param_BLP);
   
    delta_sol=output_spectral{1};
    V_sol=output_spectral{2};
```


## Detailed specification of the spectral algorithm
If we additionally specify the following parameters, we can change the settings of the spectral algorithm. If not specified, the default values are used.
For instance, if we specify spec.TOL=1e-14, the tolerance level is set to 1E-14, rather than the default value 1E-10.

* fixed_point_iter_spec (default: 1):
Let f(x) be the function specified as the input of spectral_func.
If fixed_point_iter_spec==1, we assume that f(x) is a fixed point mapping, and we solve a nonlinear equation f(x)-x=0 by the spectral algorithm.
If fixed_point_iter_spec==0, we solve for f(x)=0 by the spectral algorithm.

* TOL (default: 1e-10; scalar or vector): Tolerance level of the convergence. 
Suppose we want to solve a nonlinear equation $f(x_1,x_2)=0$.
If we use sup norm, then, we assume that an iteration converges when $||f(x_{1},x_{2})||_{\infty}<$ TOL.   
If TOL is a 1 by n_var dimensional vector, we assume that an iteration converges when $||f_{1}(x_{1},x_{2})||_{\infty}<$ TOL(1) and $||f_{1}(x_{1},x_{2})||_{\infty}<$ TOL(2) hold.


* ITER_MAX (default:1000): Maximum number of iterations

* SQUAREM_spec (default: 0)
If SQUAREM_spec==0, we apply the spectral algorithm. If SQUAREM_spec==1, alternatively apply the SQUAREM algorithm, by using SQUAREM_func.

* x_min_cell (default: []): minimum value of each type of variable. For instance, `x_min_cell={zeros(100,10),[]}` implies that the first variable we want to solve (100 by 10 dimensional array) should be nonnegative, and there is no restriction on the second variable.

* x_max_cell (default: []): maximum value of each type of variable.


* update_spec (default: []):  
If update_spec==0, apply the standard fixed point iteration. If update_spec==[], apply the spectral algorithm. Here, the step size for each type of variable is a scalar. 

* dim_hetero_alpha (default: [])
If dim_hetero_alpha>=1, we introduce variable-dimension-specific step sizes for each type of variable. For instance, suppose we want to solve a fixed point problem f(x)=x, where x is a J by T array. If we expect that the properties of the arrays `x[:,t]` (t=1,..,T) are largely different across t (time), it might be useful to introduce time-specific step sizes. If we let spec.update_spec==2, the values of x is updated by $x^{(n+1)}=x^{(n)}+\alpha^{(n)}  f(x^{(n)})$ in the spectral algorithm. Here, $\alpha^{(n)}$ is an 1 by T dimensional array. "2" implies that we introduce heterogeneity in the second dimension of the values of $\alpha^{(n)}$.
 For details, see also Section 5.3 of Fukasawa (2024).

* dampening_param (default: []):  
Suppose we want to solve f(x)=0.
In the spectral algorithm, x, the variable we want to solve, is updated by $x^{(n+1)}=x^{(n)}+\alpha^{(n)} f(x^{(n)})$.
$\alpha^{(n)}$ is computed by a prespecified formula.
If we further introduce dampening_param, the variable is alternatively updated by $x^{(n+1)}=x^{(n)}+\text{(dampening param)} \cdot \alpha^{(n)}  f(x^{(n)})$.
The introduction may lead to more stable convergence of the algorithm.


* alpha_0 (default: 1): The value of alpha_0, which should be exogenously determined in the spectral algorithm.

* alpha_max (default: 10^10): Maximum value of the step size $\alpha^{(n)}$. If the value exceeds alpha_max, let $\alpha^{(n)}=$ alpha_max.

* alpha_min (default: -10^10): Minimum value of the step size $\alpha^{(n)}$.

* common_alpha_spec (default: 0):
Suppose we want to solve a nonlinear equation $f(x_1,x_2)=0$. If common_alpha_spec==0, two types of variables x1 and x2 are updated by
$x_1^{(n+1)}=x_1^{(n)}+\alpha_{1}^{(n)} f(x_1^{(n)},x_2^{(n)})$, $x_2^{(n+1)}=x_2^{(n)}+\alpha_{2}^{(n)} f(x_1^{(n)},x_2^{(n)})$. 
If common_alpha_spec==1, we do not distinguish the type of variables, and the updating equations are $x_1^{(n+1)}=x_1^{(n)}+\alpha^{(n)} f(x_1^{(n)},x_2^{(n)})$, 
$x_2^{(n+1)}=x_2^{(n)}+\alpha^{(n)} f(x_1^{(n)},x_2^{(n)})$. See also Section 5.3 of Fukasawa (2024).
 


*  DEBUG (default: 0): If 1, display the convergence process during the iteration every 10 iterations.


* compute_alpha_spec (default: 3):  
compute_alpha_spec specifies the formula of $\alpha^{(n)}$. 
If compute_alpha_spec==1,2,3, we use the step size S1,S2,S3 defined in Fukasawa (2024). Note that they are equivalent to the three specifications in Varadhan and Roland (2008). S1 corresponds to the first step size used in Barzilai and Borwein (1988), and S2 corresponds to the second step size used in Barzilai and Borwein (1988) and the one used by La Cruz et al. (2006). Please note that some of the previous studies and built-in packages have specified the sign of the step sizes differently, and we should be careful about it. For details, see also Fukasawa (2024).
If we let compute_alpha_spec==4, we use the step size corresponding to the third specification used in the BB package (R language; Varadhan and Gilbert (2010)).  
More specifically, it is mathematically formulated as follows. Suppose we want to solve a nonlinear equation f(x)=0. Let $s^{(n)}\equiv x^{(n)}-x^{(n-1)}$, and $y^{(n)}\equiv f(x^{(n)})- f(x^{(n-1)})$.  
If compute_alpha_spec==1, let $\alpha^{(n)}=-\frac{s^{(n)\prime}y^{(n)}}{y^{(n)\prime}y^{(n)}}$.  
If compute_alpha_spec==2, let $\alpha^{(n)}=-\frac{s^{(n)\prime}s^{(n)}}{s^{(n)\prime}y^{(n)}}$.  
If compute_alpha_spec==3, let $\alpha^{(n)}=\frac{\left\Vert s^{(n)}\right\Vert _{2}}{\left\Vert y^{(n)}\right\Vert _{2}}$.  
If compute_alpha_spec==4, let $\alpha^{(n)}=sgn\left(s^{(n)\prime}y^{(n)}\right)\frac{\left\Vert s^{(n)}\right\Vert _{2}}{\left\Vert y^{(n)}\right\Vert _{2}}$.

* positive_alpha_spec (default:0)
If positive_alpha_spec==1, use positive alpha and step sizes. If $\alpha<0$, instead use 1E-8. If positive_alpha_spec==0, not restrict the sign of alpha and step sizes.

* norm_spec (default: 0): Type of norm used for assessing the convergence of the iteration.
If norm_spec==0, use sup norm. 
If norm_spec==2, use L2 norm. 
If norm_spec==10, use unit-free norm, defined by $f(x^{(n)})/x^{(n)}$, if we want to solve f(x)=0.

* Anderson_acceleration (default: 0): If 1, alternatively use (Type-II) Anderson acceleration.

## Globalization strategies (under construction)
* line_search_spec (default: 0; under construction): If 1, introduce globalization strategies, developed by La Cruz et al. (2006).
If 0, not introduce globalization strategies.  
* ITER_MAX_LINE_SEARCH (default: 10; under construction): Maximum number of iterations in each globalization step  

* M (default: 10; under construction)

* gamma (default: 0.3; under construction) 

## References
* Anderson, D. G. (1965). Iterative procedures for nonlinear integral equations. Journal of the ACM (JACM), 12(4), 547-560.
* Barzilai, J., & Borwein, J. M. (1988). Two-point step size gradient methods. IMA journal of numerical analysis, 8(1), 141-148.
* Conlon, C., & Gortmaker, J. (2020). Best practices for differentiated products demand estimation with pyblp. The RAND Journal of Economics, 51(4), 1108-1161.  
* Fukasawa, T. (2024). Fast and simple inner-loop algorithms of static/dynamic BLP estimations. arXiv preprint arXiv:2404.04494.  
* La Cruz, W., Martínez, J., & Raydan, M. (2006). Spectral residual method without gradient information for solving large-scale nonlinear systems of equations. Mathematics of computation, 75(255), 1429-1448.
* Varadhan, R., & Gilbert, P. (2010). BB: An R package for solving a large system of nonlinear equations and for optimizing a high-dimensional nonlinear objective function. Journal of statistical software, 32, 1-26.
* Varadhan, R., & Roland, C. (2008). Simple and globally convergent methods for accelerating the convergence of any EM algorithm. Scandinavian Journal of Statistics, 35(2), 335-353.


