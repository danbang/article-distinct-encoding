% Bang & Fleming (2018) Distinct encoding of decision confidence in human
% medial prefrontal cortex
%
% Visualise HDDM group-mean parameters
%
% Dan Bang 27/11/2017

% load fits
load('hddm_parameters_group_mean.mat');

% plot fits
figure('color',[1 1 1]);
ses= {'prescan','scan'};
ses_indx= [1 1 1 2 2 2];
con_indx= [4 4 4 2 2 2];
prm= {'drift','bound','nondt'};
prm_name= {'drift-rate','threshold','non-dec. time'};
prm_indx= [1 2 3 1 2 3];
ymin= [0.0 1.6 0.4 0.0 1.6 0.4];
ymax= [2.6 2.6 0.6 2.6 2.6 0.6];
xmax= [5 5 5 3 3 3];
for i_plot = 1:6
subplot(2,3,i_plot)
eval(['high_c=hddm.',ses{ses_indx(i_plot)},'.',prm{prm_indx(i_plot)},'(1:con_indx(i_plot));']);
eval(['low_c=hddm.',ses{ses_indx(i_plot)},'.',prm{prm_indx(i_plot)},'(con_indx(i_plot)+1:con_indx(i_plot)*2);']);
eval(['high_c_lb=hddm.',ses{ses_indx(i_plot)},'.lb.',prm{prm_indx(i_plot)},'(1:con_indx(i_plot));']);
eval(['low_c_lb=hddm.',ses{ses_indx(i_plot)},'.lb.',prm{prm_indx(i_plot)},'(con_indx(i_plot)+1:con_indx(i_plot)*2);']);
eval(['high_c_ub=hddm.',ses{ses_indx(i_plot)},'.ub.',prm{prm_indx(i_plot)},'(1:con_indx(i_plot));']);
eval(['low_c_ub=hddm.',ses{ses_indx(i_plot)},'.ub.',prm{prm_indx(i_plot)},'(con_indx(i_plot)+1:con_indx(i_plot)*2);']);
for i= 1:length(high_c_lb);
    plot([i-.2 i-.2],[low_c_lb(i) low_c_ub(i)],'-','color',[51/255 153/255 255/255],'LineWidth',2); hold on;
    plot([i+.2 i+.2],[high_c_lb(i) high_c_ub(i)],'-','color',[255/255 51/255 51/255],'LineWidth',2); hold on;
end
plot([1:length(low_c)]-.2, low_c,'ko','MarkerFaceColor',[51/255 153/255 255/255],'LineWidth',1); hold on
plot([1:length(high_c)]+.2, high_c,'ko','MarkerFaceColor',[255/255 51/255 51/255],'LineWidth',1); hold on
ylim([ymin(i_plot) ymax(i_plot)]); xlim([0 xmax(i_plot)]);
xlabel('distance');
ylabel('mean posterior');
set(gca,'XTiCk',1:xmax(i_plot)-1);
set(gca,'FontSize',14,'LineWidth',2);
title(prm_name{prm_indx(i_plot)},'FontWeight','normal','FontSize',16)
box('off')
end