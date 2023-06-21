function AN = construct_AN_RE(N,rho,delta,kernel)
    
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

    % % using Matlab built-in integration
    % for j=1:M+1
    %     % Interpolation with Laguerre functions
    %     AQuad(j)=integral(@(s) kernel(s).*polint(Nodes,Id(1:M+1,j),s,ww(Nodes),ones(size(s)))',-Inf,0,'AbsTol',TOL);        
    % end

    poldNDN = dNDN - rho*[zeros(M,1),eye(M)];
    polDN = DN - rho*eye(M);
    AN = polDN - (ww(Nodes(2:end))*AQuad)*(D(:,2:end)- rho*[zeros(1,M);eye(M)]); 



end