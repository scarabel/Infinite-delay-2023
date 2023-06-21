function out = PS_blowflies_inf_ext
% Copyright (c) 2023 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
%% PS_blowflies_inf_ext.m
% MatCont system definition file of the PseudoSpectral Discretization of
% the nonlinear RE defined in the paper Gurney et al, Nature, 1980,
% x(t) = beta*h(int_1^amax x(t-a)*exp(x(t-a))da)
% with h(x)= x*exp(-x)
% for the integrated state B=int_0^t x(s)ds
% The code uses the function polint.m, available from the Differentiation
% Matrix Suite by Weideman, Reddy, 2000

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= []; % @userf; % user function to select specific parameter values 
out{11}= [];
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,gamma,mu,N) 


    % construction of nodes and differentiation matrix
    rho = mu/2;
    scale = 2*rho; % for scaling of quadrature nodes and weights
    ww = @(x) exp(rho*x);
                           
    % extrema: delta = 1;
    delta = 1;
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta);

    DN = D(2:end,2:end);

    %% SYSTEM DEFINITION *** specific to the equation ***

    % Parameters and functions
    abar=1;
    beta0 = gamma*exp(mu);
    
    quad_nodes = - abar + quad_nodes;
        
%     % Without integration by parts formula
%     % compute differentiated state at quadrature nodes (for computation of right-hand side)
%     der = ww(-QuadNodes).*polint(Nodes,der_state,QuadNodes,ww(Nodes),ww(QuadNodes)); % polint function from Weideman-Reddy differentiation matrix suite
%     ARG = QuadWeights*(exp(mu*QuadNodes).*der);

    % Using integration by parts formula
    % quad_integrated = polint(Nodes,[0;state],QuadNodes,ww(Nodes),ones(size(QuadNodes)));
    quad_integrated = ww(-quad_nodes) .* polint(Nodes,[0;state],quad_nodes,ww(Nodes),ww(quad_nodes));
    ARG = quad_integrated(1)*exp(-mu) - mu * quad_weights * (exp(mu*quad_nodes).*quad_integrated);
    
%    FN = beta0*ARG.*exp(-ARG);
    FN = beta0*ARG.*exp(-100*ARG); % scaling for ease of computation: b(t)=100*btilde (equivalent for the purposes of the parameter plane)

    %% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION

    dydt= (DN - rho*eye(N))*state - ww(Nodes(2:end))*FN;

end

% --------------------------------------------------------------------------
function Weq=init(N,xeq,mu)
% INPUT: M is the discretization parameter
%        xeq is the equilibrium value of the RE
% OUTPUT Weq is the initial vector for init_EP_EP
    
    rho = mu/2;
    scale = 2*rho;
    ww = @(x) exp(rho*x);
    
    % extrema: delta = 1;
    delta = 1;
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta);

    Weq = ww(Nodes(2:end)).*Nodes(2:end)*xeq;
    
end

function out=userf(time,state,gamma,mu,N) 
% Userfunction to select specific value of parameter
    
    out = [];
    
end

% ------------
%%% AUXILIARY FUNCTIONS
