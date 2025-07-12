%%% Written by Roberto Clairmont
%%% May 2023
%%% Hikurangi PhD Research 

clear, close all; clc


%%
%load Hikurangi_Flowtest_FT.mat
load sensorpositions
load Hikurangi_Flowtest_long_FT.mat


load Sim7_Subsurface_Temp_E4_5700s_hc_1233_Detailed.mat %Generated May 16 2025 v = [8E-4 8.1E-4 8.2E-4 8.3E-4 8.4E-4 8.5E-4 8.6E-4 8.7E-4 8.8E-4 8.9E-4 9E-4 9.1E-4 9.2E-4];
load Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.mat %Generated May 16 2025 v = [8E-4 8.1E-4 8.2E-4 8.3E-4 8.4E-4 8.5E-4 8.6E-4 8.7E-4 8.8E-4 8.9E-4 9E-4 9.1E-4 9.2E-4];
load Sim7_Subsurface_Temp_E4_5700s_hc_1572_Detailed.mat %Generated May 16 2025 v = [8E-4 8.1E-4 8.2E-4 8.3E-4 8.4E-4 8.5E-4 8.6E-4 8.7E-4 8.8E-4 8.9E-4 9E-4 9.1E-4 9.2E-4];

path = 'F:\Flowmeter Analyses U1518 5 yrs'
%%

t_open_borehole = datetime('2023-03-24 18:20:00.000') %time of borehole opening

t_measured      = datetime('2023-03-24 19:55:00.000') %time of observed measurement after opening
index           = find(t == t_measured) %index of time
T_measured      = Tmap(index,:)' % Corresponding Temperature measurement


%%
% 
% totalsensors=length(FT.sn);
% t_int_FT_long  = [t(1) t(361) t(721) t(841) t(901) t(1081) t(1201) t(1411) t(1693) t(1801) t(2161) t(3241)]; %time intervals
% 
% %t_int_FT_long =(t(1):minutes(20):t(end))'; %time interval array
% T_int_FT_long = zeros(length(t_int_FT_long),totalsensors); %initalized temperature intervals array
% 
% 
%  for ii = totalsensors:-1:1
% 
%    for n = 1:length(t_int_FT_long)
% 
%      t_cell(n) = find(t == t_int_FT_long(n));
%      T_int_FT_long(n,ii) = Tmap(t_cell(n),ii);
% 
%    end
% 
%  end
% 
%  T_int_FT_long = (T_int_FT_long)'; %transpose array
%  %T_int_FT = (T_int_FT)' %transpose array
% 
%  %subplot(2,2,1)
% 
%  figure()
%  plot(t_int_FT_long,T_int_FT_long,'o','MarkerSize',10,'LineWidth',2)
%  xlabel('Time','FontSize',12)
%  ylabel('Temperature [^oC]','FontSize',12)
%  title('Temperature vs Depth at Different Time Intervals','Fontsize',16)
%  xline( [t(841) t(1693)],'-' ,{'Opening of borehole','Pull out of string'},'Fontsize',16,'linewidth',2)
%  set(gca,'FontSize', 16)
%  Screen = gcf;
%  Screen.WindowState = 'maximized';
%  pbaspect([1 1 1])
%  exportgraphics(gca,fullfile(path,'\Figures',['Temperature at Intervals.tiff']))
%  %saveas(figure(98),fullfile(path,'\Figures',['Temperature at Intervals.tiff']));

 
%% Average geotherm
load 'Tmap_2_for_avging' 
load 't_2_for_avging'
 %%

ds=60; %downsampling factor 60=10min sample period
Tmap_lowres= Tmap(1:ds:end,:);
%t_lowres= t(1:ds:end);

t_foravg = t_lowres(t_lowres >= '2023-03-01 00:05:00.000' & t_lowres <= '2023-03-23 23:55:00.000');
Tmap_lowres_foravg(:,1:30) = Tmap_lowres((t_lowres >= '2023-03-01 00:05:00.000' & t_lowres <= '2023-03-23 23:55:00.000'),1:30);

% %  Calculate statistics of Observation
T_obs_residuals      = T_measured' - Tmap_lowres_foravg;
T_obs_residuals_mean = mean(T_obs_residuals)';
T_obs_residuals_std  = std(T_obs_residuals)';

Tobs_prior_mean_per_row = (mean(Tmap_lowres_foravg))';
Tstd_prior_per_row = (std(Tmap_lowres_foravg))';

% % Calculate Best Fit line for Y_intercept and Gradient
z_model_LS = -50:1:400;
Tobs_prior_mean_per_row_range = Tobs_prior_mean_per_row([1 21],1)

FT_short = FT.z([1 21],1)

[lsm_T_int] = leastsquares(FT_short,Tobs_prior_mean_per_row_range)
g_grad = lsm_T_int(2)
yinter = lsm_T_int(1)

T_ls_model  = (lsm_T_int(1) + lsm_T_int(2).*z_model_LS)';

figure()
plot(Tobs_prior_mean_per_row,FT.z,'ob','MarkerSize',10,'MarkerFaceColor','r','DisplayName','Observation')
hold on
plot(T_ls_model,z_model_LS,'k','LineWidth',3,'DisplayName','Least Squares Fit')
title("Background geothermal gradient prior to opening of borehole to flow",'FontSize',16)
legend show Location northeast
ylabel('Depth [m]')
xlabel('Temperature [^oC]')
yline( [319.06 327.06 ],'HandleVisibility','off')
ylim([0 400])
xticks(0:2:10); 
axis square
axis ij
set(gca,'FontSize', 16)
Screen = gcf;

 %%
figure(); 
if t_measured < t_open_borehole ; % time of opening of borehole
    t_int_diff = t_open_borehole - t_measured;

    hold on
    %plot(T_int_FT_long(:,1),FT.z,'or','markersize',3,'linewidth',2)
    plot(T_measured,FT.z,'ob', 'MarkerFaceColor','g','markersize',6,'linewidth',1,'DisplayName','Observation')
    title("Temperature vs depth " + seconds(t_int_diff) + " seconds before opening",'FontSize',16)
    axis ij
    

elseif t_measured > t_open_borehole;
    t_int_diff =  t_measured - t_open_borehole;

    hold on
    plot(T_measured,FT.z,'ob','MarkerFaceColor','r','markersize',6,'linewidth',1,'DisplayName','Observation')
    plot(Tobs_prior_mean_per_row,FT.z,'ob', 'MarkerFaceColor','g','markersize',6,'linewidth',1,'DisplayName','Geotherm')

    title("Temperature vs depth " + seconds(t_int_diff) + " seconds after opening ",'FontSize',16)
    axis ij
    legend show

elseif t_measured == t_open_borehole;

    hold on
    plot(T_measured,FT.z,'ob','MarkerFaceColor','r','markersize',6,'linewidth',1,'DisplayName','Observation')
    title("Temperature vs depth during opening ",'FontSize',16)
    axis ij
end
%%

figure(97)
plot(T_measured,FT.z,'ob', 'MarkerFaceColor','r','markersize',6,'linewidth',1,'DisplayName','Observation')
title("Observed Temperature " + seconds(t_int_diff) + " seconds after opening ",'FontSize',16)
legend show Location northeast
axis ij
ylabel('Depth [m]','FontSize',12)
xlabel('Temperature [^oC]','FontSize',12)
set(gca,'FontSize', 16)
Screen = gcf;
%Screen.WindowState = 'maximized';
%exportgraphics(gca,fullfile(path,'\Figures',["Observed temperature" + seconds(t_int_diff) + " seconds after opening .tiff"]))

g_grad  = -g_grad;
T_geoth_obs = FT.z.*(-g_grad) + yinter

figure(96)
plot(T_measured - Tobs_prior_mean_per_row ,FT.z,'ob', 'MarkerFaceColor','r','markersize',6,'linewidth',1,'DisplayName','Residual: T-geoth')

legend show Location northwest
title("Observed  Residual Temperature " + seconds(t_int_diff) + " seconds after opening ",'FontSize',16)
axis ij
ylabel('Depth [m]','FontSize',12)
xlabel('Temperature [^oC]','FontSize',12)
set(gca,'FontSize', 16)
Screen = gcf;
%% Generate Geotherm for numerical model depths

inverse_Sim7_Subs_E4_Depth = 323 - Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Depth % Flip Depth values of numerical models
T_num_model_geoth = (-g_grad.*Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Depth + yinter)'; % Geotherm for numerical model of detailed depths

Residual_temp_num_models = flip(Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Temps) - T_num_model_geoth';

%% PLOT NON-RESAMPLED NUMERICAL MODEL FLOW RATES AS A FUNCTION OF DEPTH at 5700 s

v = [8E-4 8.1E-4 8.2E-4 8.3E-4 8.4E-4 8.5E-4 8.6E-4 8.7E-4 8.8E-4 8.9E-4 9E-4 9.1E-4 9.2E-4];

for n = 1:length(v)
    plotColors = turbo(size(v,2));
    ThisColor = plotColors(n,:);

figure(110)

 plot(Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Temps(:,n),inverse_Sim7_Subs_E4_Depth,'-','LineWidth',1,'Color',ThisColor,'DisplayName',""+ sprintf('%.1E' ,v(n) ) + "m/s")
 legend show Location northwest
 legend ('NumColumns',1,'FontSize',6)
 ylim([0 323])
 %pbaspect([1 1 1])
 axis square
 box on
 ylabel('Depth [m]','FontSize',12)
 xlabel('Temperature [^oC]','FontSize',12)
 title('Temperature for Mean Thermal Diffusivity [^oC]','FontSize',12')
 %hold on
 axis ij
 hold on

end
%% PLOT NON-RESAMPLED RESIDUAL NUMERICAL MODEL FLOW RATES AS A FUNCTION OF DEPTH at 5700 s

v = [8E-4 8.1E-4 8.2E-4 8.3E-4 8.4E-4 8.5E-4 8.6E-4 8.7E-4 8.8E-4 8.9E-4 9E-4 9.1E-4 9.2E-4];


for n = 1:length(v)
    plotColors = turbo(size(v,2));
    ThisColor = plotColors(n,:);

 figure(111)
 plot(Residual_temp_num_models(:,n),flip(inverse_Sim7_Subs_E4_Depth),'-','LineWidth',1,'Color',ThisColor,'DisplayName',""+ sprintf('%.1E' ,v(n) ) + "m/s")

 legend show Location northeast
 legend ('NumColumns',1,'FontSize',6)
 ylim([0 323])
 %ylim([0 400])
 %xlim([-3 4.7])
 xlim([-0.1 1])

 %pbaspect([1 1 1])
 axis square
 box on
 ylabel('Depth [m]','FontSize',12)
 xlabel('Residual Temperature [^oC] ','FontSize',12)
 title('Residual Temperature [^oC] for Mean Thermal Diffusivity [m^2/s]','FontSize',12)
 %hold on
 axis ij
 hold on

end

hold on
plot(T_measured - Tobs_prior_mean_per_row ,FT.z,'ob', 'MarkerFaceColor','r','markersize',8,'linewidth',1,'DisplayName','Residual: T-geoth')

%% RESAMPLED UNIFORM NUMERICAL MODEL DEPTHS TO NON-UNIFORM DEPTHS 

Depth_obs_0to323 = FT.z(1:21,1);

T_num_model_geoth   = (-g_grad.*Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Depth + yinter)';

T_num_model_5700s_UB    = flip(Sim7_Subsurface_Temp_E4_5700s_hc_1233_Detailed.Temps);
T_num_model_5700s_LB    = flip(Sim7_Subsurface_Temp_E4_5700s_hc_1572_Detailed.Temps);
T_num_model_5700s_Mean  = flip(Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Temps);


Resampled_T_num_model_5700s_geotherm = interp1(Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Depth,T_num_model_geoth,Depth_obs_0to323,'linear' ) ;

Resampled_T_num_model_5700s_UB      = interp1(Sim7_Subsurface_Temp_E4_5700s_hc_1233_Detailed.Depth,T_num_model_5700s_UB,Depth_obs_0to323,'linear' ) ;
Resampled_T_num_model_5700s_LB      = interp1(Sim7_Subsurface_Temp_E4_5700s_hc_1572_Detailed.Depth,T_num_model_5700s_LB,Depth_obs_0to323,'linear' );
Resampled_T_num_model_5700s_Mean    = interp1(Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Depth,T_num_model_5700s_Mean,Depth_obs_0to323,'linear' );


Resampled_temp_geotherm = interp1(Sim7_Subsurface_Temp_E4_5700s_hc_1382_Detailed.Depth,T_num_model_5700s_Mean,Depth_obs_0to323,'linear' ); 
residual_resampled_temp_models = Resampled_T_num_model_5700s_Mean - Resampled_T_num_model_5700s_geotherm;
%% PLOT RESAMPLED RESIDUAL NUMERICAL MODEL FLOW RATES AS A FUNCTION OF DEPTH at 5700 s

v = [8E-4 8.1E-4 8.2E-4 8.3E-4 8.4E-4 8.5E-4 8.6E-4 8.7E-4 8.8E-4 8.9E-4 9E-4 9.1E-4 9.2E-4];


for n = 1:length(v)
    plotColors = turbo(size(v,2));
    ThisColor = plotColors(n,:);

 figure(112)
 plot(residual_resampled_temp_models(:,n),Depth_obs_0to323,'ok','MarKerFaceColor',ThisColor,'markersize',5,'linewidth',0.5,'DisplayName',""+ sprintf('%.1E' ,v(n) ) + "m/s")

 legend show Location northeast
 legend ('NumColumns',1,'FontSize',6)
 ylim([0 323])
 xlim([-0.05 1])

 %pbaspect([1 1 1])
 axis square
 box on
 ylabel('Depth [m]','FontSize',12)
 xlabel('Residual Temperature [^oC]','FontSize',12)
 title('Resampled Numerical Model Temperatures [^oC] for Mean D_t_h [m^2/s] ')
 %hold on
 axis ij
 hold on

end

hold on
plot(T_measured - Tobs_prior_mean_per_row ,FT.z,'ob', 'MarkerFaceColor','r','markersize',8,'linewidth',1,'DisplayName','Residual_O_b_s: T-geoth')

%% REMOVAL of ROW 18 with BAD SENSOR

Resampled_T_num_model_5700s_UB(18,:) = []; 
Resampled_T_num_model_5700s_LB(18,:) = []; 
Resampled_T_num_model_5700s_Mean(18,:) = [];
Resampled_T_num_model_5700s_geotherm(18,:) = [];

T_measured = T_measured(1:21,1);
T_measured(18,:) = []; 

Tobs_prior_mean_per_row = Tobs_prior_mean_per_row(1:21,1);
Tobs_prior_mean_per_row(18,:) = [] ;

Mean_observed_residuals = T_measured - Tobs_prior_mean_per_row;


Depth_obs_0to323(18,:) = [];
Depth_obs_0to323_rem_sens_18 = Depth_obs_0to323 ;

Tstd_prior_per_row = Tstd_prior_per_row(1:21,1);
Tstd_prior_per_row(18,:) = []  ;

%% CREATE STRUCTURE WITH FLOW RATE FIELDNAMES

for i = 1:size(Resampled_T_num_model_5700s_Mean,2);

      fieldname =  sprintf('v_%g', v(i));
      fieldname = strrep(fieldname, '.','_');

     num_models_std_dev.(fieldname)    = (Resampled_T_num_model_5700s_UB(:,i) - Resampled_T_num_model_5700s_LB(:,i))/2;
     num_models_mean.(fieldname)       = Resampled_T_num_model_5700s_Mean(:,i);
     num_models_residuals.(fieldname)  = num_models_mean.(fieldname)  - Resampled_T_num_model_5700s_geotherm;

 end

%% ESTIMATING MSE FRom Monte Carlo Simultion!!!!! 05_9_2025

 num_simulations = 100000;           % Number of Monte Carlo simulations
 fit_params = zeros(num_simulations, 2);  % Store slope and intercept


 MC_obs_residuals = zeros(length(Mean_observed_residuals), num_simulations);

 for n = 1:length(Mean_observed_residuals);

     MC_obs_residuals(n,:) = Mean_observed_residuals(n) + Tstd_prior_per_row(n) .* randn(1, num_simulations);

     for ii = 1:length(fieldnames(num_models_mean));
          fieldname =  sprintf('v_%g', v(ii));
          fieldname = strrep(fieldname, '.','_');
 
          MC_num_residuals.(fieldname)(n,:) = num_models_residuals.(fieldname)(n) + num_models_std_dev.(fieldname)(n).* randn(1, num_simulations);
          MC_residual_differences_loop.(fieldname)(n,:) = MC_num_residuals.(fieldname)(n,:) - MC_obs_residuals(n,:);
          square_residual_differences_loop.(fieldname)(n,:) = MC_residual_differences_loop.(fieldname)(n,:).^2;
          mse_loop.(fieldname) = mean(square_residual_differences_loop.(fieldname));       
     end
    
     allMSEs = struct2cell(mse_loop);
     allMSEs = cell2mat(allMSEs);
     [MinValue, MinIndex] = min(allMSEs);

     % Map the index to the flow rate vector
     Minimized_flowrate = v(MinIndex);

 end

 Mean_flowrate = mean(Minimized_flowrate)
 Std_Dev_flowrate = std(Minimized_flowrate)
 Confidence_Int_flowrate = prctile(Minimized_flowrate, [2.5, 97.5])
 %% Histogram

 figure();
 %histogram(Minimized_flowrate,13,'EdgeColor','black')
 histogram(Minimized_flowrate,'BinWidth',0.0000094)
 xlim([min(v)  max(v)])
 ticklocations = linspace(min(v),max(v),length(v)); 
 ticks = ticklocations(1:2:end)
 xticks(ticks)
 pbaspect([1 1 1])

 % yticklocations = linspace(1,length(Minimized_flowrate),length(Minimized_flowrate));
 % ticks = yticklocations(1:9999:end)
 % yticks(ticks)
 xlabel('Flow Rate (m/s)')
 ylabel('Number of Realizations')
 axis square
 set(gca,'FontSize', 16)
 exportgraphics(gca,fullfile(path,'\Figures',["Histogram of MC Simulation for FlowRates.pdf"]))



 % end


