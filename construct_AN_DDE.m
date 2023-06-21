function AN = construct_AN_DDE(N,rho,delta,kernel,a)
    
    ww = @(x) exp(rho.*x);       
    Id = eye(N+1);
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta);

    M = length(Nodes)-1;
    DN = D(2:end,2:end);
    dNDN = D(2:end,:);
    
    AQuad=zeros(1,M+1); 
    % using quadrature (Quarteroni, Sacco, Saleri)
    for j=1:M+1
        % Interpolation with Laguerre functions
        AQuad(j)= quad_weights*(kernel(quad_nodes).*polint(Nodes,Id(1:M+1,j),quad_nodes,ww(Nodes),ones(size(quad_nodes))));
    end
    AQuad(1) = AQuad(1) + a;
    poldNDN = dNDN - rho*[zeros(M,1),eye(M)];
    polDN = DN - rho*eye(M);
    AN = [AQuad; poldNDN];
end