%% linear_tests.m
% copyright Francesca Scarabel 2023
% Code to reproduce the tests of the linear iDEs for the paper
% Scarabel, Vermiglio, "Equations with infinite delay:\\ pseudospectral
% discretization for numerical stability and bifurcation in an abstract framework"
% The following codes use the suites by:
% Gautschi, Gauss–Radau formulae for Jacobi and Laguerre weight functions, Math. Com- put. Simulation, 54 (2000), pp. 403–412, 1999 International Symposium on Computational Sciences, to honor John R. Rice (West Lafayette, IN).
% Weideman, Reddy, A MATLAB differentiation matrix suite, ACM T.Math.Software, 26 (2000), 465–519.

%% Test on linear DDE
% script to approximate the eigenvalues of the linear DDE with infinite
% delay:
% y'(t) = a y(t) + int_(-Inf)^0 k(s) y(t+s) ds
% and to study the convergence behavior

clear; 
close all; 

savefigure = 0; % set to 1 to save figures
colormap = lines(10);

% Definitions of parameters for each test
% Uncomment the desidered set of parameters

% Test a1 (exponential, lambda=0)
parset = 'a1';
mu = 2;
a = 3; 
k0 = -a*mu; 
kernel = @(theta) k0*exp(mu*theta);
Lambda = 0;
% 
% % Test a2 (exponential, same as a1 but lambda = (a-mu)
% parset = 'a2';
% mu = 2;
% a = 3; 
% k0 = -a*mu; 
% kernel = @(theta) k0*exp(mu*theta);
% Lambda = a-mu;
% 
% % Test b (exponential, imaginary, lambda = \pm i*mu)
% parset = 'b';
% mu = 2;
% a = 2; 
% k0 = -2*mu.^2; 
% kernel = @(theta) k0*exp(mu*theta);
% Lambda = 1i*mu; % char roots computed by Laplace transform

% % Test c (exponential, real)
% parset = 'c';
% mu = 2;
% a = 3*mu; 
% k0 = -(a+mu)^2/4; 
% kernel = @(theta) k0*exp(mu*theta);
% Lambda = (a-mu)/2; % char roots computed by Laplace transform

% % Test d (gamma)
% parset = 'd';
% a = 0; 
% k0 = []; 
% shape_gamma = 2;
% rate_gamma = 4; % mean_gamma = shape_gamma/rate_gamma;
% mu = rate_gamma;
% kernel = @(theta) gampdf(-theta,shape_gamma,1/rate_gamma);
% Lambda = fsolve(@(x) (x-a).*(x+mu).^2-mu^2, 0); % char roots computed by Laplace transform

% % Test e (gamma)
% parset = 'e';
% % parset = 'e_int'; TOL = 1e-12; % comment this line if the integral is computed via quadrature
% a = 0; 
% k0 = []; 
% shape_gamma = 4.5; % pi; % exp(1);
% rate_gamma = 4; % mean_gamma = shape_gamma/rate_gamma;
% mu = rate_gamma;
% kernel = @(theta) gampdf(-theta,shape_gamma,1/rate_gamma);
% Lambda = fsolve(@(x) x-a- integral(@(y) kernel(-y).*(exp(-x*y)),0,Inf), 1) % char roots computed by Laplace transform
% % Lambda = [];


% Definition of parameters for discretization
% Uncomment the desired set of parameters

rho = mu/2; rho_label='rho';
% rho = mu/4; rho_label='rho2';
nodes_type = 'zeros'; delta = 0;
% nodes_type = 'ext'; delta = 1;

% Test name
nametest = ['LinDDE_',nodes_type,'_',rho_label,'_',parset];

% Routine
Nmax = 100; % 100;
Nlist=[[1:9],[10:2:19],[20:5:49],[50:10:Nmax]]; % Nlist for plotting purposes

% For-loop to compute the errors varying N
for N = Nmax:-1:1
    
    % display(['Calculating index ',num2str(iN)])    
    % Computation of standard Laguerre nodes
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho,delta);

    AN = construct_AN_DDE(N,rho,delta,kernel,a);
    [EigenVectorAN,LambdaAN] = eig(AN);
    LambdaAN = diag(LambdaAN);
    
    % The following if statement is only used for using the built-in
    % integral of Matlab
    if strcmp(parset,'e_int')
        AN = construct_AN_DDE_integral(N,rho,delta,kernel,a,TOL);
        [EigenVectorAN,LambdaAN] = eig(AN);
        LambdaAN = diag(LambdaAN);
    end

    if isempty(Lambda) && N==Nmax
        [~,i] = max(real(LambdaAN));
        Lambda = LambdaAN(i);
    end
    [error_eig,ind] = min(abs(Lambda-LambdaAN));
    LambdaApprox(N,:) = LambdaAN(ind);
    Error_Eigs(N) = error_eig;
    
    % Approximation of exponential eigenfunctions
    weigfcn = @(x) exp((rho+Lambda)*x); % reference eigenfunction
    wMaxError(N) = max(abs(weigfcn(Nodes)-EigenVectorAN(:,ind)./EigenVectorAN(1,ind))); % error compared to normalised eigenvector (s.t. equal 1 in theta=0)

end

fhDDE = figure; clf;
loglog(Nlist,Error_Eigs(Nlist),'.','LineStyle','-','Color','k','HandleVisibility','on'); hold on %colormap(testtype,:)
loglog(Nlist,wMaxError(Nlist),'.','LineStyle',':','Color','k','HandleVisibility','off'); hold on
grid on
xlabel('$N$','interpreter','latex'); 
% ylabel('error eig','interpreter','latex');
axis([1 100 1e-16 1])
title(strcat('$a=$',num2str(a),', $k_0= $',num2str(k0),', $\mu= $',num2str(mu),', $\lambda=$ ',num2str(Lambda,4)),'interpreter','latex');
legend([nodes_type,' $\rho=$',num2str(rho)], 'interpreter','latex')

% The next if statement is to estimate the convergence order for the
% polynomial convergence
if strcmp(parset,'e')
    conv_order = log(Error_Eigs(100)/Error_Eigs(50))/(log(100/50))
    % loglog(Nlist,Nlist.^conv_order,'-r','HandleVisibility','off')
end

if savefigure
    savefig(fhDDE,[pwd '/Figures/',nametest,'.fig']);
    saveas(fhDDE,[pwd '/Figures/',nametest],'png');
    display('Figures saved');
end

%% Test on linear RE

% script to approximate the eigenvalues of the linear RE with infinite
% delay:
% y(t) = int_(-Inf)^0 k(s) y(t+s) ds
% and to study the convergence behavior

close all; 
clear;

savefigure = 0; % set to 1 to save figures
colormap = lines(10);

% Definitions of parameters for each test
% Uncomment the desired set of parameters

% Test f (Sinusoidal)
parset = 'h';
mu = 1;
a = 1;
k0 = 1; 
kernel = @(theta) k0*(sin(-a*theta)+1).*exp(mu*theta);
Lambda = []; % COMPUTE REFERENCE VALUE
Pol = [1;
    3*mu-k0;
    3*mu^2-2*k0*mu+a^2-k0*a;
    mu^3-k0*mu^2+mu*a^2-k0*a^2-k0*a*mu];
Lambdavec = roots(Pol)
Lambda = Lambdavec(end); 

% % Test g (Sinusoidal)
% parset = 'i';
% mu = 1.5;
% a = 1;
% k0 = 3; % mu*(mu^2+a^2)/(mu^2+a^2+a*mu); 
% kernel = @(theta) k0*(sin(-a*theta)+1).*exp(mu*theta);
% Lambda = []; % COMPUTE REFERENCE VALUE
% Pol = [1;
%     3*mu-k0;
%     3*mu^2-2*k0*mu+a^2-k0*a;
%     mu^3-k0*mu^2+mu*a^2-k0*a^2-k0*a*mu];
% Lambdavec = roots(Pol)
% Lambda = Lambdavec(1); 

% % Test h (exponential, lambda=0)
% parset = 'f';
% mu = 2;
% k0 = mu; 
% kernel = @(theta) k0*exp(mu*theta);
% Lambda = k0-mu;

% % Test i (exponential, lambda=k0-mu)
% parset = 'g';
% mu = 2;
% k0 = 5; 
% kernel = @(theta) k0*exp(mu*theta);
% Lambda = k0-mu;

% Definition of parameters for discretization (select two combinations)
% Uncomment the desired set of parameters

rho1 = mu/2; rho_label='rho';
% rho1 = mu/3; rho_label='rho2';
nodes_type = 'zeros'; delta = 0;
% nodes_type = 'ext'; delta = 1;

rho = mu; % for weighted L1 norm (w.r.t rho)

% Test name
nametest = ['LinRE_',nodes_type,'_',rho_label,'_',parset];

% Routine
Nmax = 100;
Nlist=[[1:9],[10:2:19],[20:5:49],[50:10:Nmax]]; % Nlist for plotting purposes

% For-loop to compute the errors varying N
for N = Nmax:-1:1
    
    % display(['Calculating index ',num2str(iN)])    
    % Computation of standard Laguerre nodes
    [Nodes,D,quad_nodes,quad_weights] = PSD_laguerre_standard_nodes(N,rho1,delta);

    AN = construct_AN_RE(N,rho1,delta,kernel);
    [EigenVectorAN,LambdaAN] = eig(AN);
    LambdaAN = diag(LambdaAN);

    if isempty(Lambda) && N==Nmax
        [~,i] = max(real(LambdaAN));
        Lambda = LambdaAN(i);
    end    

    [error_eig,ind] = min(abs(Lambda-LambdaAN));
    LambdaApprox(N,:) = LambdaAN(ind);
    Error_Eigs(N) = error_eig;
    
    rho1_EigenVectorAN_diff = (D(:,2:end)- rho1*[zeros(1,N);eye(N)])*EigenVectorAN(:,ind); % consider original state (derivative of eigenvector)
    rho1_EigenVectorAN_diff = rho1_EigenVectorAN_diff./rho1_EigenVectorAN_diff(1); % normalisation such that is 1 in theta=0.
    rho_EigenVectorAN_diff = exp((rho-rho1)*Nodes).*rho1_EigenVectorAN_diff; % consider original state (derivative of eigenvector)

    % Approximation of exponential eigenfunctions
    weigfcn = @(x) exp((rho+Lambda)*x); % weighted reference eigenfunction

    wL1Error(N) = norm(abs(weigfcn(Nodes)-rho_EigenVectorAN_diff),1); % discrete 1-norm of vector
    
end

fhRE = figure; clf;
loglog(Nlist,Error_Eigs(Nlist),'.','LineStyle','-','Color','k','HandleVisibility','on'); hold on %colormap(testtype,:)
% loglog(Nlist,wMaxError(Nlist),'.','LineStyle','--','Color','k','HandleVisibility','off'); hold on
loglog(Nlist,wL1Error(Nlist),'.','LineStyle',':','Color','k','HandleVisibility','off'); hold on
grid on
xlabel('$N$','interpreter','latex'); 
% ylabel('error eig','interpreter','latex');
axis([1 100 1e-16 1])
title(strcat('$k_0= $',num2str(k0),', $\mu= $',num2str(mu),', $\lambda=$ ',num2str(Lambda,4)),'interpreter','latex');
legend([nodes_type,' $\rho=$',num2str(rho1)], 'interpreter','latex')

if savefigure
    savefig(fhRE,[pwd '/Figures/',nametest,'.fig']);
    saveas(fhRE,[pwd '/Figures/',nametest],'png');
    display('Figures saved');
end

