function out = PS_beretta_breda_ext
% Copyright (c) 2023 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 

out{1} = @init;
out{2} = @fun_eval;
out{3} = []; %@jacobian;
out{4} = []; %@jacobianp;
out{5} = []; %@hessians;
out{6} = []; %@hessiansp;
out{7} = []; %@der3;
out{8} = [];
out{9} = [];
out{10}= []; %@userf; % user function to select specific parameter values 
out{11}= [];
out{12}= [];
out{13}= [];
end

% --------------------------------------------------------------------------
function dydt = fun_eval(time,state,ttau,deathA,deathJ,a,b,m,N) 

    % construction of nodes and differentiation matrix
    
    rho = 0.25*(deathJ+m/ttau);
    ww = @(x) exp(rho*x);
    
    % extrema: delta = 1;
    delta = 1;
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta);

    DN = D(2:end,2:end);
    dNDN = D(2:end,:);

    %% SYSTEM DEFINITION *** specific to the equation ***
                                   
    % Parameters and functions
    Gamma = gampdf(-quad_nodes,m,ttau/m);
%    Gamma =
%    (n/ttau)^n.*(-quad_nodes).^(n-1)/prod(1:n-1).*exp(quad_nodes*n/ttau); %
%    don't use this - built-in function is better
       
%    state_unscaled = subplus(state([3:N]))./ww(quad_nodes); % setting to zero negative values is particularly important for infinite delay
    % When using Radau-Laguerre quadrature:
    state_unscaled = subplus(state)./ww(quad_nodes); % setting to zero negative values is particularly important for infinite delay
    
    % For more general quadrature nodes, may need interpolation:
    %    state_unscaled1 = subplus(polint(Nodes,state,quad_nodes,ww(Nodes),ww(quad_nodes)))./ww(quad_nodes);

    Kernel = exp(deathJ*quad_nodes-a*state_unscaled);
    ARG = Gamma.*Kernel.*state_unscaled;
   
    GM = - deathA*state(1) + b*(quad_weights*ARG);

    %% FINAL APPROXIMATING ODE SYSTEM - PSEUDOSPECTRAL DISCRETIZATION

    dydt= [GM; (dNDN-rho*[zeros(N,1),eye(N)])*state];

end
 
 
% --------------------------------------------------------------------------
function Weq=init(N,yeq,ttau,deathA,deathJ,a,b,m) 
% INPUT: M is the discretization parameter
%        xeq,yeq are column vectors with the equilibrium states of the RE and
%        DDE respectively
% OUTPUT Weq is the initial vector for init_EP_EP
    
    rho = 0.25*(deathJ+m/ttau);
    ww = @(x) exp(rho*x);
    
    % extrema: delta = 1;
    delta = 1;
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta);

    Weq=ww(Nodes)*yeq;
    
end


% function out=userf(time,state,gamma,mu,aux,tau,M) 
% % Userfunction to select specific value of parameter
%     
%     out = (gamma-60).*(gamma-70);
%     
% end

% ------------
%%% AUXILIARY FUNCTIONS
