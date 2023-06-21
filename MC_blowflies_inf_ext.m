% Copyright (c) 2023 Francesca Scarabel
% This code is distributed under the MIT license, see LICENSE for 
% licensing information. 
% 
% MC_blowflies_inf_ext
% command line instructions for the MatCont continuation of the system defined
% in PS_blowflies_inf_ext.m
% (Gurney et al, Nature 1980)

%%

clear;
clearvars -global cds
close all

savefigure = 0;
continue_periodic = 0;

%% Analytic Hopf curve

w = [1.7:0.01:3]'; % linspace(pi/2,pi-0.1,100)';

c1= w.*cos(w)./sin(w);
c2= - w./sin(w);

mmu = - c1;
bbeta = - c1.*exp(1+c2./c1);

figure(2)
plot(mmu,bbeta./mmu,'k','LineWidth',2); hold on
plot([0;5], [1;1], 'g'); hold on % existence
axis([0 5 0 50])
xlabel('mu'); ylabel('beta/mu')

%% Numerical continuation with Matcont

% Discretization parameters
% N=20;
% Continuation parameters
ap1=1; % index of the continuation parameter in the vector par
ap2=2;
ap3=3;
TOL=1e-6;
TestTOL=1e-6;

for N = [10, 20, 30] % iN = [1,2,3]
    
    % N = 10*iN;
    display(['Running index ',num2str(N)])


% Initial parameter values
mu=3;
gamma=1;
%gamma = 10;
par=[gamma,mu,N]';

% Approximated equilibrium corresponding to par
xeq=0; 
%xeq = 2;

%% Continuation process

MM=N; % dimension of the approximating ODE system
handles=feval(@PS_blowflies_inf_ext);
opt=contset;
global cds;

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

state_eq=feval(handles{1},N,xeq,mu); % initializes equilibrium vector
[x0,v0]=init_EP_EP(@PS_blowflies_inf_ext,state_eq,par0,ap1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=0;
while ((length(se)<3) && xe(end,end)>0.1 && jj<10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
    jj=jj+1;
end

figure(1); clf(figure(1));
xlabel('$\gamma$','interpreter','latex','fontsize',12);
ylabel('state','interpreter','latex','fontsize',12);
cpl(xe,ve,se,[MM+1 1]);

%% Detection of singular points
% xe,ve,se,he,fe = output of previous continuation
% par = current parameter vector
% ap1 = index of continuation parameter in vector par

% BP, branching point
for ii=1:length(se)
    if strcmp(se(ii).label,'BP')==1
        BP_index=se(ii).index;
        sBP=se(ii);
        break;
    end
end
par(ap1)=xe(end,BP_index);
BP=xe(1:MM,BP_index);

xeBP=xe; veBP=ve; seBP=se; heBP=he; feBP=fe;
parBP=par;

%% Equilibrium continuation from BP
% BP = vector of variables at BP
% parBP = parameter vector at BP
% sBP = information about BP from previous continuation
display('Starting equilibrium continuation from BP');

% set options
%opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'MaxStepsize',1);
opt=contset(opt,'Backward',0);

[x0,v0]=init_BP_EP(@PS_blowflies_inf_ext,BP,parBP,sBP,0.1);
[xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
jj=0;
while ( xe(end,end)<30 && xe(end,end)>0 && jj<10)
    [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
    jj=jj+1;
end

% Plot
figure(1) 
cpl(xe,ve,se,[MM+1 1]);
xlabel('$\gamma$','interpreter','latex','fontsize',12);
title(['Bifurcation Example 2.2, mu=',num2str(mu),', M=',num2str(N)]);

%% Detection of singular points
% H, Hopf point

for ii=1:size(se)
    if (strcmp(se(ii).label,'H ')==1 && xe(end,se(ii).index)>0)
        H_index=se(ii).index;
        break;
    end
end
par(ap1)=xe(end,H_index);
H=xe(1:MM,H_index);

xeH=xe; veH=ve; seH=se; heH=he; feH=fe;
parH=par;

%% H continuation in two parameters
% H = vector of variables at H
% parH = parameter vector at H
% ap1,ap2 = index of continuation parameters in the vector par
display('Starting H continuation');

% set options
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxStepsize',0.1);

% forward branch
opt=contset(opt,'Backward',0);

[x0,v0]=init_H_H(@PS_blowflies_inf_ext,H,parH,[ap1 ap2]);
[xh,vh,sh,hh,fh]=cont(@hopf,x0,[],opt); xh(MM+1,end)

jj=0;
while (xh(MM+2,1)<5 && xh(MM+2,end)<5 && xh(MM+2,1)>0.6 && xh(MM+2,end)>0.6 && jj<5)
     [xh,vh,sh,hh,fh]=cont(xh,vh,sh,hh,fh,cds); xh(MM+1,end)
      jj=jj+1;
end

% % For M = 10
% par(ap1)=xh(MM+1,sh(2).index);
% par(ap2)=xh(MM+2,sh(2).index);
% H2=xh(1:MM,sh(2).index);
% parH2=par;
% 
% opt=contset(opt,'Backward',0);
% [x0,v0]=init_H_H(@PS_blowflies_inf_ext,H2,parH2,[ap1 ap2]);
% [xh2,vh2,sh2,hh2,fh2]=cont(@hopf,x0,[],opt); xh2(MM+1,end)


%% backward branch
opt=contset(opt,'Backward',1);

[x0,v0]=init_H_H(@PS_blowflies_inf_ext,H,parH,[ap1 ap2]);
[xhb,vhb,shb,hhb,fhb]=cont(@hopf,x0,[],opt); xhb(MM+1,end)
% [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)

jj=0;
while (xhb(MM+2,1)<5 && xhb(MM+2,end)<5 && xhb(MM+2,1)>0.6 && xhb(MM+2,end)>0.6 && jj<5)
    [xhb,vhb,shb,hhb,fhb]=cont(xhb,vhb,shb,hhb,fhb,cds); xhb(MM+1,end)
    jj=jj+1;
end

% Plot stability regions
figure(2); hold on
fig = plot(xh(MM+2,:),xh(MM+1,:)./xh(MM+2,:),'HandleVisibility','off'); % Hopf
col1 = get(fig,'color');
plot(xhb(MM+2,:),xhb(MM+1,:)./xhb(MM+2,:), 'Color',col1)
axis([0 5 0 50]);

figure(3); hold on
plot(xh(MM+2,:),xh(MM+1,:)./xh(MM+2,:),'b','HandleVisibility','off'); % Hopf
plot(xhb(MM+2,:),xhb(MM+1,:)./xhb(MM+2,:),'b')
xlabel('mu'); ylabel('gamma/mu');
title(['Stability regions, M=',num2str(N)]);
% savefig(['M',num2str(M),'_PD_curve_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]); % ]); %
if savefigure
    savefig(gcf,[pwd '/Figures/BF_ext_Hopf_N_',num2str(N),'.fig'])
    saveas(gcf,[pwd '/Figures/BF_ext_Hopf_N_',num2str(N)],'png')
%    matlab2tikz([pwd '/Figures/BF_ext_Hopf_N_',num2str(N)]);
    display('Figures saved');
end

% Plot of bifurcation diagram
figure(1); clf; % bifurcation diagram in the sun-star variable
cpl(xeBP,veBP,seBP,[MM+1 1]); hold on;
cpl(xeH,veH,seH,[MM+1 1]);
xlabel('gamma');
title(['Bifurcation Example 2.2, mu=',num2str(mu),', N=',num2str(N)]);
% savefig(['M',num2str(M),'_bif_diagram_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);


if continue_periodic
%% Limit cycle continuation from H
% H = vector of variables at H
% parH = parameter vector at H
% ap1 = index of continuation parameter in vector par
display('Starting LC continuation from H');

% These settings do not work for M = 50
% set options
TOL=1e-6; 
opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL);

opt=contset(opt,'MaxNumPoints',200);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'Multipliers',1);
opt=contset(opt,'MaxStepsize',1); %100
opt=contset(opt,'InitStepsize',1e-1);
opt=contset(opt,'Adapt',1); % opt=contset(opt,'Adapt',0); % for M=10

ntst=40; % number of intervals
ncol=4; % degree of polynomial

[x0,v0]=init_H_LC(@PS_blowflies_inf_ext,H,parH,ap1,1e-1,ntst,ncol);
[xlc,vlc,slc,hlc,flc]= cont(@limitcycle,x0,v0,opt); xlc(end,end)
jj=1;
while ((length(slc)<5) && xlc(end,end)< 90 && jj<10)
    [xlc,vlc,slc,hlc,flc]= cont(xlc,vlc,slc,hlc,flc,cds); xlc(end,end)
    jj=jj+1;
end

% save(['M',num2str(M),'_Example22_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);

%% detection of PD bifurcation

SPD=1;
for S=1:size(slc)
    if strcmp(slc(S).label,'PD ')==1
        PD_index=slc(S).index;
        SPD=S;
        break;
    end
end

par_PD=xlc(end,PD_index);
parPD=parH;
parPD(ap1)=par_PD;

%% PD continuation in two parameters % ap1,ap2 = index of continuation
% parameters in the vector par display('Starting PD continuation');
% 
TOL=1e-6;
FunTOL=1e-6;

% set options
opt=contset(opt,'FunTolerance',FunTOL); opt=contset(opt,'VarTolerance',TOL);
opt=contset(opt,'TestTolerance',TOL); %opt=contset(opt,'Singularities',1);
opt=contset(opt,'Eigenvalues',0);
opt=contset(opt,'Multipliers',1); 
opt=contset(opt,'Singularities',0);
opt=contset(opt,'MaxNumPoints',100);
opt=contset(opt,'MaxStepsize',1);

% forward branch
opt=contset(opt,'Backward',0);

[x0,v0]=init_PD_PD(@PS_blowflies_inf_ext,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
[xpd,vpd,spd,hpd,fpd]=cont(@perioddoubling,x0,v0,opt); % 

    [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); 
    [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); 

l=size(xpd,1);

% jj=0; 
% while (xpd(l,1)<5 && xpd(l,end)<5 && xpd(l,1)>0.6 && xpd(l,end)>0.6 && jj<5) %
%     [xpd,vpd,spd,hpd,fpd]=cont(xpd,vpd,spd,hpd,fpd,cds); 
%     xpd(l,end) %   
%     xpd(end-1,end)./xpd(end,end)
%     jj=jj+1;
% end

% save(['M',num2str(M),'_Example22_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);

%% backward branch
opt=contset(opt,'Backward',1);
opt=contset(opt,'MaxStepsize',1);
% opt=contset(opt,'Adapt',3);

[x0,v0]=init_PD_PD(@PS_blowflies_inf_ext,xlc,slc(SPD),[ap1 ap2],ntst,ncol);
[xpd1,vpd1,spd1,hpd1,fpd1]=cont(@perioddoubling,x0,v0,opt); % 
    [xpd1,vpd1,spd1,hpd1,fpd1]=cont(xpd1,vpd1,spd1,hpd1,fpd1,cds); xpd1(l,end) %       
    [xpd1,vpd1,spd1,hpd1,fpd1]=cont(xpd1,vpd1,spd1,hpd1,fpd1,cds); xpd1(l,end) %       
% 
% jj=0; 
% while (xpd1(l,1)<4.1 && xpd1(l,end)<4.1 && xpd1(l,1)>0.6 && xpd1(l,end)>0.6 && jj<5) %  
%     [xpd1,vpd1,spd1,hpd1,fpd1]=cont(xpd1,vpd1,spd1,hpd1,fpd1,cds); xpd1(l,end) %       
%     jj=jj+1;
% end

% save(['M',num2str(M),'_Example22_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]);

who('global')
global MC cds eds hds lds sOutput
save(['M',num2str(N),'_blowflies_ntst_',num2str(ntst)]);

figure(2); hold on
plot(xpd(end,:),xpd(end-1,:)./xpd(end,:),'r'); % Period doubling
plot(xpd1(end,:),xpd1(end-1,:)./xpd1(end,:),'r'); % Period doubling

figure(1);
plot(xlc(end,:),upperbound,'g',xlc(end,:),lowerbound,'g');
% ylabel('max/min','interpreter','latex')

for ii=2:length(slc)-1
    index=slc(ii).index;
    plot(xlc(end,index),upperbound(index),'or',xlc(end,index),lowerbound(index),'or');
end

end


end

figure(2); 
xlabel('mu'); ylabel('gamma/mu');
legend('reference','existence','$N=10$','$N=20$','$N=30$','Interpreter','latex')
% savefig(['M',num2str(M),'_PD_curve_tau_',num2str(tau_max),'_ntst_',num2str(ntst)]); % ]); %
if savefigure
    savefig(gcf,[pwd '/Figures/BF_ext_plane.fig'])
    saveas(gcf,[pwd '/Figures/BF_ext_plane.png'],'png')
%    matlab2tikz([pwd '/Figures/BF_ext_plane.tex']);
    display('Figures saved');
end


%% Error analysis on Hopf

close all
%%

TOL=1e-10;
TestTOL=1e-10;

ap1=1; % index of the continuation parameter in the vector par

Nmax = 51;
HopfN = nan(Nmax,1);
BranchPointN = nan(Nmax,1);
LambdaBPN = nan(Nmax,1);
LambdaHN = nan(Nmax,1);

% Initial parameter values
gamma=1;
xeq=0; 
handles=feval(@PS_blowflies_inf_ext);

% parameter and reference values
mu=2; 
H1 = []; %24.852286382216239; % calculated with complete interpolation and N=50
% 24.852614133380420; N=30
% Hopf(end);  
BP1 = 2; % BranchPoint(end);

% % parameter and reference values
% mu=3; 
% H1 = []; % 29.692815712722695; % calculated with complete interpolation and N=50
% % Hopf(end);  
% BP1 = 3; % BranchPoint(end);

% % parameter and reference values
% mu=4; 
% H1 = []; % 35.693476826645309; % calculated with complete interpolation and N=50
% % Hopf(end);  
% BP1 = 4; % BranchPoint(end);

LambdaBP = [];
LambdaH = [];

for N = 1:Nmax
    
    display(['Running index ',num2str(N)])
    
    par=[gamma,mu,N]';
    MM=N; % dimension of the approximating ODE system
    
    opt=contset;
    global cds;

    % Equilibrium continuation from initial point [xeq;yeq]
    display('Starting equilibrium continuation from initial point');
    par0=par;

    % set options
    opt=contset(opt,'Singularities',1);
    opt=contset(opt,'FunTolerance',TOL); opt=contset(opt,'VarTolerance',TOL);
    opt=contset(opt,'TestTolerance',TestTOL);
    opt=contset(opt,'Eigenvalues',1);
    opt=contset(opt,'Backward',0);
    opt=contset(opt,'MaxNumPoints',200);

    state_eq=feval(handles{1},N,xeq,mu); % initializes equilibrium vector
    [x0,v0]=init_EP_EP(@PS_blowflies_inf_ext,state_eq,par0,ap1);
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    jj=0;
    while ((length(se)<3) && xe(end,end)>0.1 && jj<10)
        [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds);
        jj=jj+1;
    end

    figure(1); clf(figure(1));
    xlabel('$\gamma$','interpreter','latex','fontsize',12);
    ylabel('state','interpreter','latex','fontsize',12);
    cpl(xe,ve,se,[MM+1 1]);

    % BP, branching point
    for ii=1:length(se)
        if strcmp(se(ii).label,'BP')==1
            BP_index=se(ii).index;
            sBP=se(ii);
            break;
        end
    end
    par(ap1)=xe(end,BP_index);
    BP=xe(1:MM,BP_index);

    xeBP=xe; veBP=ve; seBP=se; heBP=he; feBP=fe;
    parBP=par;

    BranchPointN(N) = xe(end,BP_index);
    LambdaBPN(N) = fe(end,BP_index);
    
    % Equilibrium continuation from BP
    % BP = vector of variables at BP
    % parBP = parameter vector at BP
    % sBP = information about BP from previous continuation
    display('Starting equilibrium continuation from BP');

    % set options
    %opt=contset(opt,'MaxNumPoints',100);
    opt=contset(opt,'MaxStepsize',0.1);
    opt=contset(opt,'Backward',0);

    [x0,v0]=init_BP_EP(@PS_blowflies_inf_ext,BP,parBP,sBP,0.1);
    [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
    jj=0;
    while ( xe(end,end)<30 && xe(end,end)>0 && jj<10)
        [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
        jj=jj+1;
    end

    % Plot
    figure(1) 
    cpl(xe,ve,se,[MM+1 1]);
    xlabel('$\gamma$','interpreter','latex','fontsize',12);
    title(['Bifurcation Example 2.2, mu=',num2str(mu),', M=',num2str(N)]);

    % Detection of singular points
    % H, Hopf point

    for ii=1:size(se)
        if (strcmp(se(ii).label,'H ')==1 && xe(end,se(ii).index)>0)
            if max(real(fe(:,se(ii).index-1)))<0 && max(real(fe(:,se(ii).index+1)))>0 
                HopfN(N) = xe(end,se(ii).index);
                LambdaHN(N) = fe(end,se(ii).index);
                break;
            end
        end
    end
    
    if ~isfinite(HopfN(N))
        opt=contset(opt,'MaxStepsize',0.05);
        [x0,v0]=init_BP_EP(@PS_blowflies_inf_ext,BP,parBP,sBP,0.1);
        [xe,ve,se,he,fe]=cont(@equilibrium,x0,v0,opt); xe(end,end)
        jj=0;
        while ( xe(end,end)<40 && xe(end,end)>0 && jj<10)
            [xe,ve,se,he,fe]=cont(xe,ve,se,he,fe,cds); xe(end,end)
            jj=jj+1;
        end
        for ii=1:size(se)
            if (strcmp(se(ii).label,'H ')==1 && xe(end,se(ii).index)>0)
                if max(real(fe(:,se(ii).index-1)))<0 && max(real(fe(:,se(ii).index+1)))>0 
                    HopfN(N) = xe(end,se(ii).index);
                    break;
                end
            end
        end
    end
    
end

% Plot error
% global MC cds eds sOutput  
% save(['Error_Hopf_TOL10_mu',num2str(mu)])

% Using as reference value N=50:
if isempty(H1)
    BP1 = BranchPointN(end);
    H1 = HopfN(end);
    LambdaBP = LambdaBPN(end);
    LambdaH = LambdaHN(end);
end

figure;
loglog(1:Nmax,abs(BranchPointN-BP1),'.-','LineWidth',1,'Color','k'); hold on
loglog(1:Nmax,abs(HopfN-H1),'o-','LineWidth',1,'Color','k'); hold on
loglog(1:Nmax,abs(LambdaBPN-LambdaBP),'.:','LineWidth',1,'Color','k','HandleVisibility','on')
loglog(1:Nmax,abs(LambdaHN-LambdaH),'o:','LineWidth',1,'Color','k','HandleVisibility','on')
legend('BP','H','lambda BP','lambda H')
xlabel('$N$','Interpreter','latex')
xlim([0 Nmax])
grid on
title(['$H=$',num2str(H1),],'Interpreter','latex')

if savefigure
    save(['BF_ext_Error_Hopf_TOL10_mu',num2str(mu),'.mat'],'H1','BP1','LambdaBP','LambdaH','mu','gamma','-mat');
    savefig(gcf,[pwd '/Figures/BF_ext_Error_Hopf_TOL10_mu',num2str(mu),'.fig'])
    saveas(gcf,[pwd '/Figures/BF_ext_Error_Hopf_TOL10_mu',num2str(mu)],'png')
%    matlab2tikz([pwd '/Figures/BF_ext_Error_Hopf_TOL10_mu',num2str(mu)]);
    display('Figures saved');
end