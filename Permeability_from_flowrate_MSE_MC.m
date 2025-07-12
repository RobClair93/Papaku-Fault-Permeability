clc; close all; clear
path = 'F:\Flowmeter Analyses U1518 5 yrs'
del_P   = [3.05E6];
time_perm = 0:1:864000; %10 days
k_perm = [4E-15 4.5E-15 5E-15 5.5E-15 6.0E-15 6.5E-15 7.E-15 8E-15 8.5E-15 9E-15...
          1E-14 1.1E-14 1.2E-14 1.3E-14 1.4E-14 1.5E-14 1.6E-14 1.7E-14 1.8E-14 1.9E-14 2E-14];
%k_perm = [1E-15 1.65E-14 1E-14 1E-13 1E-12 1E-11];
phi     = 0.45;
H       = 8;
Beta    = 5E-10;
th_diff = 3.28E-07;
dyn_vis = 1.38*10^-3;
r       = 0.165% radius of screen interval %0.0514;
gamma   = 0.572
for nn = 1:length(k_perm)
tao_perm(:,nn) = (k_perm(nn).*time_perm)./(phi.*dyn_vis.*Beta.*r^2);


end

for i = 1:length(k_perm) %orig

        for ii = 1:length(tao_perm) %orign
  
    
    vflow(ii,i)  =  ( (8.*k_perm(i).*H.*del_P)./(pi.^2.*dyn_vis.*r.^2) )* ( 1./(log(4.*tao_perm(ii,i))-(2.*gamma))  -  gamma./(abs(log(4.*tao_perm(ii,i))-(2.*gamma))).^2     ) ;


                    %eqn ater Fisher et al 1997                    
        end 
figure(101)
loglog(time_perm,vflow(:,i),'DisplayName',"k = " + k_perm(i) + "m^2 ",'LineWidth',3) %orig
legend show Location southwest
title("Fluid velocity vs Time for Excess Pressure \DeltaP =" + del_P/1E6 +"MPa",'FontSize',16) %orig
xlabel('log_1_0Time [s]','FontSize',16)
ylabel('Fluid velocity [m/s]','FontSize',16)
xlim([25 10E5])
axis square
hold on
% loglog(5700,1.2E-3,'ob','MarkerFaceColor','g','markersize',8,'HandleVisibility','off')
% loglog(554400,1E-4,'ob','MarkerFaceColor','g','markersize',8,'HandleVisibility','off')
loglog(5700,8.6E-4,'ob','MarkerFaceColor','g','markersize',8,'HandleVisibility','off')
loglog(5700,8.81E-4,'ob','MarkerFaceColor','y','markersize',8,'HandleVisibility','off')
loglog(5700,8.39E-4,'ob','MarkerFaceColor','r','markersize',8,'HandleVisibility','off')

loglog(554400,1.83E-4,'ob','MarkerFaceColor','g','markersize',8,'HandleVisibility','off')
loglog(554400,2.4E-4,'ob','MarkerFaceColor','y','markersize',8,'HandleVisibility','off')
loglog(554400,2.97E-4,'ob','MarkerFaceColor','r','markersize',8,'HandleVisibility','off')

xline((5700),'-',{'5700 s'},'Fontsize',12,'linewidth',2,'HandleVisibility','off')
xline((554400),'-',{'6.42 d'},'Fontsize',12,'linewidth',2,'HandleVisibility','off')
set(gca,'fontsize',16)

set(findall(gcf,'type','text'),'fontsize',16)

end
%%


vobs1 = [8.6E-4 - 0.21E-4, 8.6E-4 ,8.6E-4 + 0.21E-4];

vflow_all_at_5700 = zeros(size(vobs1,2),size(k_perm,2));
vflow_all_at_5700_diff = zeros(size(vobs1,2),size(k_perm,2));
vflow_all_at_5700_square  = zeros(size(vobs1,2),size(k_perm,2));
vflow_all_at_5700_min     = zeros(size(vobs1,2),1);

for i= 1:length(vobs1);

vflow_all_at_5700(i,:) = vflow(5700,:);

vflow_all_at_5700_diff(i,:)   = vobs1(i) - vflow_all_at_5700(i,:);
vflow_all_at_5700_square(i,:) = vflow_all_at_5700_diff(i,:).^2;
vflow_all_at_5700_min(i)  = min(vflow_all_at_5700_square(i,:));

[MinValue(i), MinIndex(i)] = min(vflow_all_at_5700_square(i,:));


end

min_perm1 = k_perm(MinIndex)
%%

vobs2 = [2.4E-4 - 0.56E-4, 2.4E-4 ,2.4E-4 + 0.56E-4];

vflow_all_at_154hrs = zeros(size(vobs2,2),size(k_perm,2));
vflow_all_at_154hrs_diff = zeros(size(vobs2,2),size(k_perm,2));
vflow_all_at_154hrs_square  = zeros(size(vobs2,2),size(k_perm,2));
vflow_all_at_154hrs_min     = zeros(size(vobs2,2),1);

for i= 1:length(vobs2);

vflow_all_at_154hrs(i,:) = vflow(5700,:);

vflow_all_at_154hrs_diff(i,:)   = vobs2(i) - vflow_all_at_154hrs(i,:);
vflow_all_at_154hrs_square(i,:) = vflow_all_at_154hrs_diff(i,:).^2;
vflow_all_at_154hrs_min(i)  = min(vflow_all_at_154hrs_square(i,:));

[MinValue(i), MinIndex(i)] = min(vflow_all_at_154hrs_square(i,:));


end

min_perm2 = k_perm(MinIndex)


