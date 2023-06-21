function [coll_nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta)
% INPUT:
% N: discretization index
% rho: scaling parameter
% delta: either Laguerre zeros (delta=0) or extrema (delta=1)
% nodes
% OUTPUT:
% coll_nodes: N+1 collocation nodes
% D: full differentiation matrix D-\rho I
% quad_nodes, quad_weights: quadrature rules, delta=0 corresponds to GL formula (N-2
% nodes); delta=1 corresponds to GRL formula (N-1 nodes including zero)

    unscaled_nodes = zeros(N+1,1);
    scale = 2*rho;
   
    % Computation of standard Laguerre nodes
    if delta==0
        ab0=r_laguerre(N,delta);
        xw = gauss(N,ab0);
        unscaled_nodes(2:end) = xw(:,1);
        quad_nodes = - xw(:,1)/(2*rho);
        quad_weights = (xw(:,2).*exp(xw(:,1))/(2*rho))'; % weights for unweighted integration 
        % alternative for zeros--laguerre:
        % zeros_laguerre = lagroots(N-2); % using Weideman-Reddy differentiation suite
    elseif delta==1
        % GRL rules with additional nodes 0
        xw = radau_laguerre(N,0);
        unscaled_nodes(2:end) = xw(2:end,1);
        quad_nodes = - xw(:,1)/(2*rho);
        quad_weights = (xw(:,2).*exp(xw(:,1))/(2*rho))'; % weights for unweighted integration 
    end
    D_unscaled = poldif(unscaled_nodes, exp(-0.5*unscaled_nodes), -0.5*ones(1,N+1));

    coll_nodes = - unscaled_nodes/scale;
    D = - scale*D_unscaled;

end