function [Mmodel] = tectono_magmatic_fcn(zm,pm_rate,vpl,PARS)

% a simple toy model of tectono-magmatic interactions at a mid-ocean ridge
% J.-A. Olive & P. Dublanchet, July 2019 ? July 2020
% all units SI

dsxxf=PARS.mu*(PARS.rho*(1-PARS.la)*PARS.g*zm)/(sqrt(1+PARS.mu^2)+PARS.mu); %--fault slip event threshold
dsxxdm=PARS.T+PARS.psf+PARS.rho*PARS.g*zm;        %--max eruption threshold (assumes minimal magma pressure is 0, when all magma has drained from AML)
dsige=PARS.E*PARS.w_dyke/PARS.w_rift;  %--eruption stress drop
dsigs=PARS.E*PARS.slip/PARS.w_rift;   %--earthquake stress drop

niter=5e5;     %--nb iterations
dt=0.01*dsxxdm/pm_rate;  %--time step

t=zeros(1,niter);
sxxt=zeros(1,niter);sxxt(1)=0;
sxxd=zeros(1,niter);sxxd(1)=dsxxdm;
te=[];
ts=[];
flg=0;

iter=1;
while iter<=niter

    sxxd(iter+1)=sxxd(iter)-pm_rate*dt; %--pressure builds up, diking threshold decreases
    sxxt(iter+1)=sxxt(iter)+dt*(PARS.E/PARS.w_rift)*vpl; %--axial tension builds up due to far-field plate separation
    t(iter+1)=t(iter)+dt;
    
    
    if sxxd(iter+1)<dsxxf %--stress is below tectonic threshold
        if sxxt(iter+1)>=sxxd(iter+1) %--dike intrusion occurs
           dtcorr=(sxxd(iter)-sxxt(iter))*dt/(sxxt(iter+1)-sxxt(iter)-sxxd(iter+1)+sxxd(iter));
           t(iter+1)=dtcorr+t(iter);
           sxxd(iter+1)=sxxd(iter)-pm_rate*dtcorr;
           sxxt(iter+1)=sxxt(iter)+dtcorr*(PARS.E/PARS.w_rift)*vpl;
           
           iter=iter+1;
           t(iter+1)=t(iter);

           sxxd(iter+1)=dsxxdm; %--magma pressure resets to zero, dike threshold becomes maximal
           sxxt(iter+1)=max([0 sxxt(iter)-dsige]); %--axial stress drops by amount dsige (intrusion stress drop)
           te=[te t(iter+1)];
           
        end
    else
        if sxxt(iter+1)>=dsxxf %--fault slip event occurs
           dtcorr=(dsxxf-sxxt(iter))*dt/(sxxt(iter+1)-sxxt(iter));
           t(iter+1)=dtcorr+t(iter);
           sxxd(iter+1)=sxxd(iter)-pm_rate*dtcorr;
           sxxt(iter+1)=sxxt(iter)+dtcorr*(PARS.E/PARS.w_rift)*vpl;
           
           iter=iter+1;
           t(iter+1)=t(iter);
           
           sxxd(iter+1)=sxxd(iter); % diking threshold not affected by earthquake
           sxxt(iter+1)=max([0 sxxt(iter)-dsigs]); %--axial stress drops by amount dsigs (earthquake / slow slip stress drop)
           ts=[ts t(iter+1)];

        end
    end
    
    if sxxt(iter+1)==0 && flg==0 % stress has dropped below zero (goes compressive), which is not possible 
        % this can occur when far-field stresses build up too quickly (e.g., when w_rift is too small)
        disp('ERROR - AXIAL STRESS DROPPED TO ZERO')
        disp('SUM OF EQs AND DIKES WILL NOT MATCH FAR-FIELD')
        flg = 1;
    end
    
    iter=iter+1;
end

if isempty(te) || isempty(ts)
    Mmodel = NaN;
else
    ue=PARS.w_dyke*[1:1:length(te)]; % total amount of horizontal displacement due to dike
    us=PARS.slip*[1:1:length(ts)]; % total amount of horizontal displacement due to fault slip
  
    us=interp1(ts,us,t);
    ue=interp1(te,ue,t);
   
    M = ue./(ue+us);
    Mmodel = min(M(isnan(M)==0)); 
end

if flg==1
    Mmodel = NaN;
end


plot_cycle=1;
if plot_cycle==1
% plot stress cycle
figure();clf
plot(t/(365.25*24*3600)-5e3,sxxt/1e6,'-k','linewidth',1.2)
hold on
plot(t/(365.25*24*3600)-5e3,sxxd/1e6,'-','Color',[0.7 0 0],'linewidth',1.2);
plot(t/(365.25*24*3600)-5e3,dsxxf*ones(size(t))/1e6,'--k','linewidth',1.2);
xlabel('time (years)','fontsize',15);
ylabel('stress (MPa)','fontsize',15);
ylim([0 220])
xlim([0 500])
grid on
legend('axial stress','diking threshold','tectonic threshold')
set(gca,'PlotBoxAspectRatio',[2 1 1])
end

