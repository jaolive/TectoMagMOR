clear 
% a simple analytical model of tectono-magmatic interactions at a mid-ocean ridge
% predicting M vs. spreading rate
% J.-A. Olive & P. Dublanchet, July 2019 ? July 2020
% all units SI

%% LITHOSPHERE THICKNESS 
% 2 end-member scenarios for lithospheric thickness vs. spreading rate
regime = 'thin'; % 'thin', for thin lithosphere endmember; or 'thick' for thick lithosphere endmember

%% WHICH MODEL TO USE
yn_sim=0; % set to 0 to use the analytical model (quick); 1 to use the toy model, can take a while unless N_rate and N_press are reduced to ~ 100 and 100.


%% MODEL PARAMETERS
yr=365.25*24*3600;  %--year to sec
rho=2.9e3;          %--density of rock (kg/m3)
rhow=1000;          %--density of water (kg/m3)
la=rhow/rho;        %--hydrostatic pressure ratio in the crust around the fault
g=9.81;             %--gravity (m/s^2)
psf=1000*g*2000;    %-- pressure at the seafloor (Pa)
w_rift=1e4;         %--axis width (m)
w_dyke=1;           %--dyke width (m)
slip=1;             %--earthquake coseismic slip (m)
mu=0.6;             %--friction coeff
E=10e9;             %--Young's modulus (Pa)
nu=0.25;            %--Poisson's ratio
T=1e6;              %--Crustal tensile strength (Pa)
Hcru=6000;          %--Crustal thickness (m)
G = E/(2*(1+nu));   %--Shear modulus (Pa)

%% SCENARIOS FOR PRESSURE RATE VS. SPREADING RATE
% extreme values of ?Pmdot/?U = ximax and ximin
ximax = (4*Hcru./(pi*(1-nu)))*(12e9/(10*50^2));
ximin = (4*Hcru./(pi*(1-nu)))*(2e9/(100*500^2));


%% DEFINE GRID OF PARAMETERS 
N_rate=250; % number of spreading rate values
N_press=250; % number of pressure rate values




%% CONSTRUCT GRIDS OF PRESSURE RATE, SPREADING RATE, LITHOSPHERE THICKNESS 
pprate = linspace(.005,.4,N_press);
vvpl = linspace(.0005,0.18,N_rate)/24/3600/365.25;
[vpl,pm_rate]= meshgrid(vvpl,pprate);
% a model relating AMC depth to spreading rate based on observations
if strcmp(regime,'thin')==1
    amcdepth = @(vpl) 1200+(6000-1200)*exp(-(vpl*100*24*365.25*3600)/1.8); % works well AND acts as envelope for AMC data
    col='r';
elseif strcmp(regime,'thick')==1
    %amcdepth = @(vpl) 1200+(22000-1200)*exp(-(vpl*100*24*365.25*3600)/3); % to have deeper AML at ultraslow
    amcdepth = @(vpl) 2000+(26000-2000)*exp(-(vpl*100*24*365.25*3600)/2.7); % to have deeper AML at ultraslow
    col='b';
end
zm=amcdepth(vpl);
prate_model = sqrt(ximax*ximin)*vvpl; % intermediate scenario for pmdot vs. spreading rate

%% STORE PARAMETERS IN ARRAY
PARS.yr = yr;
PARS.rho = rho;
PARS.rhow = rhow;
PARS.la = la;
PARS.g = g;
PARS.psf = psf;
PARS.w_rift = w_rift;
PARS.w_dyke = w_dyke;
PARS.slip = slip;
PARS.mu = mu;
PARS.E = E;
PARS.T = T;

 
%% ANALYTICAL MODEL
dsxxf=mu*(rho*(1-la)*g*zm)/(sqrt(1+mu^2)+mu); %--fault slip event threshold
dsxxdm=T+psf+rho*g*zm;        %--max eruption threshold (assumes minimal magma pressure is 0, when all magma has drained from AML)
dsige=E*w_dyke/w_rift;  %--eruption stress drop
dsigs=E*slip/w_rift;   %--earthquake stress drop
tau_er = (dsxxdm-dsxxf)./pm_rate; %-- eruption recurrence time
M = (w_dyke./(tau_er.*vpl)); % M value
M(M<0)=0;
M(M>1)=1;
tau_s  = slip./(vpl.*(1-M)); % earthquake recurrence time


%% SIMULATIONS (if using the toy model)
if yn_sim==1
Msim = zeros(size(zm));
cnt=0;
for i = 1:length(pprate)
    for j = 1:length(vvpl)
        % full model of tectono-magmatic cycles
        [Msim(i,j)] = tectono_magmatic_fcn(zm(i,j),pm_rate(i,j),vpl(i,j),PARS); 
        cnt=cnt+1;
        if mod(cnt,10)==0
            disp(cnt)
        end
    end
end
Msim(isnan(Msim)==1)=1;
end


figure(3) % CONTOUR PLOT OF M VS. PRESSURE RATE AND SPREADING RATE
hh=pcolor(100*vvpl*365.25*24*3600,pprate,M);

if yn_sim==1
    hold on
    epsil = 1e-4;
    contour(100*vvpl*365.25*24*3600,pprate,Msim,[epsil .1 .2 .3 .4 .5 .6 .7 .8 .9 1-epsil],'k')
end

set(hh,'Edgecolor','none')
xlabel('spreading rate (cm/yr)');
ylabel('magma pressure rate (Pa/s)');
set(gca,'PlotBoxAspectRatio',[1 1 1])
colorbar
caxis([0 1])
M_model = interp2(vpl,pm_rate,M,vvpl,prate_model);
hold on
plot(100*vvpl*365.25*24*3600,prate_model,'k')
plot(100*vvpl*365.25*24*3600,ximax*vvpl,'k--')
plot(100*vvpl*365.25*24*3600,ximin*vvpl,'k--')



% DATA
[Md, ML, MU, U, mor] = textread('M_data.txt','%f%f%f%f%s'); % M vs. SPREADING RATE
[AMCdepth, AMCdepthrange, U_AMC, name_AMC]=textread('depth_of_AMLs.txt','%f%f%f%s'); % AML DEPTH vs. SPREADING RATE
[EQmin, EQmax, U_EQ]=textread('depth_of_microEQs.txt','%f%f%f'); % MICRO-EQ DEPTH vs. SPREADING RATE
%load seismic_crustal_thickness





%% PAPER FIGURES

figure(1) % M DATASET
I=find(strcmp(mor,'MAR')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'bo')
grid on
hold on
I=find(strcmp(mor,'EPR')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'rp')
I=find(strcmp(mor,'GSC')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'bx')
I=find(strcmp(mor,'ELSC')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'kd')
I=find(strcmp(mor,'JDF')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'c^')
I=find(strcmp(mor,'SEIR')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'cs')
I=find(strcmp(mor,'CR')==1);
errorbar(U(I),Md(I),Md(I)-ML(I),MU(I)-Md(I),'cv')
ylabel('M');
xlabel('spreading rate (cm/yr)');
set(gca,'PlotBoxAspectRatio',[1.5 1 1])
xlim([0 18])
ylim([.3 1])


figure(2) % THICKNESS, PRESSURE RATE, AND M ACROSS RATES
subplot(311)
errorbar(U_AMC,1e-3*AMCdepth,1e-3*AMCdepthrange/2,'ko')
grid on
hold on
plot(U_EQ,EQmax,'ks')
axis ij
plot(100*vvpl*365.25*24*3600,amcdepth(vvpl)/1e3,col)
ylabel('axial thickness (km)');
set(gca,'PlotBoxAspectRatio',[1.5 1 1])
xlim([0 18])
ylim([0 17])

subplot(312)
plot(vvpl*365.25*24*3600*100,prate_model,'k')
hold on
plot(vvpl*365.25*24*3600*100,ximax*vvpl,'k--')
plot(vvpl*365.25*24*3600*100,ximin*vvpl,'k--')
plot([11 11],[1.9e-1 4.6e3],'k','linewidth',2)
plot([11 11],[4e-5 2.3e-3],'k','linewidth',2)
grid on
set(gca,'PlotBoxAspectRatio',[1.5 1 1])
xlim([0 18])
ylim([1e-5 1e4])
set(gca,'Yscale','log')
ylabel('pressure build-up rate (Pa/s)');
xlabel('spreading rate (cm/yr)');

subplot(313)
errorbar(U,Md,Md-ML,MU-Md,'ko')
grid on
hold on
plot(100*vvpl*365.25*24*3600,M_model,col)
ylabel('M');
xlabel('spreading rate (cm/yr)');
set(gca,'PlotBoxAspectRatio',[1.5 1 1])
xlim([0 18])



%% transition thickness
mup = mu/(mu+sqrt(1+mu^2));
gam = 1-mup*(1-la);
P0 = T+psf;



thick = linspace(.1,20e3,1000);
ximed = sqrt(ximin*ximax);
xiM = sqrt(sqrt(ximin*ximax)*sqrt(ximed*ximax));
xim = sqrt(sqrt(ximin*ximax)*sqrt(ximed*ximin));
MvalueMED = w_dyke*sqrt(ximin*ximax)./(P0+gam*rho*g*thick);
MvalueMIN = w_dyke*ximin./(P0+gam*rho*g*thick);
MvalueMAX = w_dyke*ximax./(P0+gam*rho*g*thick);
MvalueMEDplus = w_dyke*sqrt(ximed*ximax)./(P0+gam*rho*g*thick);
MvalueMEDminus = w_dyke*sqrt(ximed*ximin)./(P0+gam*rho*g*thick);
MM = w_dyke*xiM./(P0+gam*rho*g*thick);
Mm = w_dyke*xim./(P0+gam*rho*g*thick);
xivec = logspace(log10(ximin),log10(ximax),1000);

Hcrit = (w_dyke*xivec-P0)./(gam*rho*g);
HcritDetach = (w_dyke*xivec/0.5-P0)./(gam*rho*g);



figure(5) % EFFECT OF LITHOSPHERE THICKNESS ON M, THICKNESS THRESHOLDS
subplot(121)
plot(thick/1e3,min(MvalueMED,1),'k','linewidth',2)
xlabel('lithospheric thickness (km)');
ylabel('M');
set(gca,'PlotBoxAspectRatio',[1 1 1])
grid on
ylim([0 1])
hold on
plot(thick/1e3,min(MvalueMAX,1),'k','linewidth',1)
plot(thick/1e3,min(MvalueMIN,1),'k','linewidth',1)
plot(thick/1e3,min(MvalueMEDplus,1),'k--','linewidth',1)
plot(thick/1e3,min(MvalueMEDminus,1),'k--','linewidth',1)
plot(thick/1e3,min(MM,1),'k--','linewidth',1)
plot(thick/1e3,min(Mm,1),'k--','linewidth',1)

subplot(122)
plot(xivec,Hcrit/1e3,'k','linewidth',2)
hold on
plot(xivec,HcritDetach/1e3,'b','linewidth',2)
set(gca,'PlotBoxAspectRatio',[1 1 1])
grid on
hold on
set(gca,'Xscale','log')
ylim([0 10])
ylabel('transition thickness (km)');
xlabel('\chi');
plot([sqrt(ximin*ximax) sqrt(ximin*ximax)],[0 10],'k--','linewidth',2)
plot([sqrt(ximin*ximed) sqrt(ximin*ximed)],[0 10],'k--')
plot([sqrt(ximed*ximax) sqrt(ximed*ximax)],[0 10],'k--')
plot([xiM xiM],[0 10],'k--')
plot([xim xim],[0 10],'k--')
plot([ximin ximin],[0 10],'k--')
plot([ximax ximax],[0 10],'k--')


show_regimes=1;
if show_regimes==1

%% PLOT SEISMO-VOLCANIC CYCLES IN FOUR EXAMPLE SITUATIONS
[Ss] = tectono_magmatic_fcn(3000,sqrt(ximin*ximax)*(0.025/365.25/24/3600),0.025/365.25/24/3600,PARS);
ylim([16 108])
xlim([0 250])
title('slow, symmetric')

[Sa] = tectono_magmatic_fcn(8000,sqrt(ximin*ximax)*(0.025/365.25/24/3600),0.025/365.25/24/3600,PARS);
ylim([47 250])
xlim([0 250])
title('slow, asymmetric')

[I] = tectono_magmatic_fcn(2500,sqrt(ximin*ximax)*(0.06/365.25/24/3600),0.06/365.25/24/3600,PARS);
ylim([10 93])
xlim([0 250])
title('intermediate')

[F] = tectono_magmatic_fcn(1700,sqrt(ximin*ximax)*(0.15/365.25/24/3600),0.15/365.25/24/3600,PARS);
ylim([0 70])
xlim([0 250])
title('fast')

end



