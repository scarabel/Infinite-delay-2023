% Copyright (c) 2023 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 

%%
clear;
clearvars -global cds
close all

%% Numerical continuation with Matcont

savefigure = 0;

TOL=1e-6;
TestTOL=1e-6;

for N = [5, 10, 20] % iN = [1,2,3]
    
    % N = 10*iN;
    display(['Running index ',num2str(N)])
    NN = N+1; % dimension of the approximating ODE system
    % Initial parameter values
    ttau=9;
    deathA=0.5; deathJ=1; a=7; b=350; m=7;
    par=[ttau, deathA, deathJ, a, b, m, N]'; % parameter vector

    % Approximated equilibrium corresponding to par
    yeq=0.1; % tau=9;

    % Continuation parameters
    ap1=1; % index of the bifurcation parameter in the vector 'par'
    ap2=6;

    %% Continuation process
    handles=feval(@PS_beretta_breda_ext);
    opt=contset;
    global cds;

    state_eq=feval(handles{1},N,yeq,ttau,deathA,deathJ,a,b,m); % initializes equilibrium vector

    % %%
    % % 
    % rhs = @(t,y) feval(handles{2},t,y,ttau,deathA,deathJ,a,b,n,M);
    % 
    % [tt,yy] = ode45(rhs,[0 100],state_eq);
    % plot(tt,yy(:,1))
    % 
    % % state_eq = yy(end,:)

    %% Equilibrium continuation from initial point [xeq;yeq]

    display('Starting equilibrium continuation from initial point');
    par0=par;

    % set options
    opt=contset(opt,'Singularities',1);
    opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
    opt=contset(opt,'TestTolerance',TestTOL);
    opt=contset(opt,'Eigenvalues',1);
    opt=contset(opt,'Backward',0);
    opt=contset(opt,'MaxNumPoints',200);

    [x0,v0]=init_EP_EP(@PS_beretta_breda_ext,state_eq,par,ap1);
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    % jj=0;
    % while ((length(se)<3) && xe(end,end)>0.1 && jj<10)
    %     [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
    %     jj=jj+1;
    % end


    figure(1); clf(figure(1));
    xlabel('$\tau$','interpreter','latex','fontsize',12);
    ylabel('$\overline{N}$','interpreter','latex','fontsize',12);
    cpl(xe,ve,se,[NN+1 1]);

    %% Detection of singular points
    % H, Hopf point

    for ii=1:size(se)
        if (strcmp(se(ii).label,'H ')==1 && xe(end,se(ii).index)>0)
            H_index=se(ii).index;
            break;
        end
    end

    par(ap1)=xe(end,H_index);
    H=xe(1:NN,H_index);

    xeH=xe; veH=ve; seH=se; heH=he; feH=fe;
    parH=par;

    %% H continuation in two parameters
    % H = vector of variables at H
    % parH = parameter vector at H
    % ap1,ap2 = index of continuation parameters in the vector par
    display('Starting H continuation');

    % set options
    % TOL = 1e-6;
    opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
    opt=contset(opt,'TestTolerance',TOL);
    opt=contset(opt,'Singularities',0);
    % opt=contset(opt,'InitStepsize',1e-6);
    % opt=contset(opt,'MaxStepsize',1e4);
    opt=contset(opt,'MaxNumPoints',200);
    % opt=contset(opt,'Adapt',1);

    % forward branch
    opt=contset(opt,'Backward',0);

    [x0,v0]=init_H_H(@PS_beretta_breda_ext,H,parH,[ap1 ap2]);
    [xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(NN+1,end)
    % [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+2,end)

    % backward branch
    opt=contset(opt,'Backward',1);

    [x0,v0]=init_H_H(@PS_beretta_breda_ext,H,parH,[ap1 ap2]);
    [xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(NN+1,end)
    % [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)

%     figure(2); hold on
%     cpl(xh,vh,sh,[MM+1 MM+2]);
%     cpl(xhb,vhb,shb,[MM+1 MM+2]);
    
    figure(2)
    fig = plot(xh(NN+1,:),xh(NN+2,:)); hold on
    c1 = get(fig,'color');
    plot(xhb(NN+1,:),xhb(NN+2,:),'Color',c1,'HandleVisibility','off')

    % save(['Hopf_curve_N',num2str(N)]);
end

axis([ 1 4 0 12])
xlabel('$\tau$','Interpreter','latex')
ylabel('$n$','Interpreter','latex')
legend('$N=5$','$N=10$','$N=20$','Interpreter','latex')

if savefigure
    savefig(gcf,[pwd '/Figures/BB_plane_ext.fig']);
    saveas(gcf,[pwd '/Figures/BB_plane_ext'],'png');
%    matlab2tikz([pwd '/Figures/BB_plane_ext']);
    display('Figures saved');
end

%% Error analysis on Hopf

% close all

deathA=0.5; deathJ=1; a=7; b=350;

% Reference values for Hopf bifurcations obtained with the Matcont
% continuation of the ODE equivalent system (file MC_ODE_beretta_breda.m)
% using tolerance TOL = 1e-10

% case 1: n=7
m = 7;
H1 = 3.130292408159863;
H2 = 1.551923477991856;
Lambda1 = 0.000000000022915 - 0.973435637011011i;
Lambda2 = -0.000000000133379 - 1.523341372437741i;

% case 2: n=7.5
m = 6.5;
H1 = [];
H2 = [];
Lambda1 = [];
Lambda2 = [];

TOL=1e-10;
TestTOL=1e-10;

Nmax = 31;
Hopf1N = nan(Nmax,1);
Hopf2N = nan(Nmax,1);
Lambda1N = nan(Nmax,1);
Lambda2N = nan(Nmax,1);
Time_Eval = nan(Nmax,1);
Time_Cont = nan(Nmax,1);

for N = 1:Nmax
    
    display(['Running index ',num2str(N)])
    NN=N+1; % dimension of the approximating ODE system
    % Initial parameter values
    ttau=9; yeq=0.1; %tau=9;
    par=[ttau, deathA, deathJ, a, b, m, N]'; % parameter vector

    % Continuation parameters
    ap1=1; % index of the bifurcation parameter in the vector 'par'
    ap2=6;

    %% Continuation process
    handles=feval(@PS_beretta_breda_ext);
    state_eq=feval(handles{1},N,yeq,ttau,deathA,deathJ,a,b,m); % initializes equilibrium vector

    % rhs = @(t,y) feval(handles{2},t,y,ttau,deathA,deathJ,a,b,n,M);
    % [tt,yy] = ode45(rhs,[0 100],state_eq);
    % plot(tt,yy(:,1))
    % state_eq = yy(end,:)
    
    opt=contset;
    global cds;

    %% Timing of one right-hand-side evaluation
    tic
    feval(handles{2},[],state_eq,ttau,deathA,deathJ,a,b,m,N);
    Time_Eval(N) = toc;
    
    %% Equilibrium continuation from initial point [xeq;yeq]

    display('Starting equilibrium continuation from initial point');
    par0=par;

    % set options
    opt=contset(opt,'Singularities',1);
    opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
    opt=contset(opt,'TestTolerance',TestTOL);
    opt=contset(opt,'Eigenvalues',1);
    opt=contset(opt,'Backward',0);
    opt=contset(opt,'MaxNumPoints',200);

    [x0,v0]=init_EP_EP(@PS_beretta_breda_ext,state_eq,par,ap1);
    tic
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    Time_Cont(N) = toc;

    figure(1); clf(figure(1));
    xlabel('$\tau$','interpreter','latex','fontsize',12);
    ylabel('$\overline{N}$','interpreter','latex','fontsize',12);
    cpl(xe,ve,se,[NN+1 1]);

    if (strcmp(se(2).label,'H ')==1)
        Hopf1N(N)= xe(end,se(2).index);
        Lambda1N(N) = fe(end,se(2).index);
        if strcmp(se(3).label,'H ')==1
            Hopf2N(N)= xe(end,se(3).index);
            Lambda2N(N) = fe(end,se(3).index);
        end
    
    else
        opt=contset(opt,'Backward',1);
        [x0,v0]=init_EP_EP(@PS_beretta_breda_ext,state_eq,par,ap1);
        [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
        if (strcmp(se(2).label,'H ')==1)
            Hopf1N(N)= xe(end,se(2).index);
            Lambda1N(N) = fe(end,se(2).index);
            if strcmp(se(3).label,'H ')==1
                Hopf2N(N)= xe(end,se(3).index);
                Lambda2N(N) = fe(end,se(3).index);
            end
        end
    end
    
end

% Using as reference value N=50:
if isempty(H1)
    H1 = Hopf1N(end);
    H2 = Hopf2N(end);
    Lambda1 = Lambda1N(end);
    Lambda2 = Lambda2N(end);
end

% Plot error
loglog(1:Nmax,abs(Hopf1N-H1),'.-','LineWidth',1,'Color','k'); hold on
loglog(1:Nmax,abs(Hopf2N-H2),'o-','LineWidth',1,'Color','k','HandleVisibility','on')
loglog(1:Nmax,abs(Lambda1N-Lambda1),'.:','LineWidth',1,'Color','k','HandleVisibility','on')
loglog(1:Nmax,abs(Lambda2N-Lambda2),'o:','LineWidth',1,'Color','k','HandleVisibility','on')
grid on
xlabel('$N$','Interpreter','latex')
xlim([0 Nmax])
legend('H1','H2','$\lambda_1$','$\lambda_2$','Interpreter','latex')
title('Error on Hopf','Interpreter','latex')


if savefigure
    savefig(gcf,[pwd '/Figures/BB_ext_Error_Hopf_ext_TOL10_m',num2str(m),'.fig']);
    saveas(gcf,[pwd '/Figures/BB_ext_Error_Hopf_ext_TOL10_m',num2str(m)],'png');
%    matlab2tikz([pwd '/Figures/BB_ext_Error_Hopf_ext_TOL10_m',num2str(m)]);
    display('Figures saved');
end

figure
plot(1:Nmax,Time_Cont,'.-','LineWidth',1,'Color','k'); hold on
% plot(1:Nmax,Time_Eval,'o-','LineWidth',1,'Color','k'); hold on
grid on
xlabel('$N$','Interpreter','latex')
ylabel('seconds','Interpreter','latex')
xlim([0 Nmax])
legend('complete GRL','Interpreter','latex')
title('Time for a 100-point continuation','Interpreter','latex')

if savefigure
    savefig(gcf,[pwd '/Figures/BB_ext_Time_ext_TOL10_m',num2str(m),'.fig']);
    saveas(gcf,[pwd '/Figures/BB_ext_Time_ext_TOL10_m',num2str(m)],'png');
%     matlab2tikz([pwd '/Figures/BB_ext_Time_ext_TOL10_m',num2str(m)]);
    display('Figures saved');
end