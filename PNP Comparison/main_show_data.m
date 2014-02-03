clear
close all
% % load data_P3P_C3
% load data_QuanLan_C3
% % save data_QuanLan_C3 diff_T_epnp diff_Rot_epnp diff_T_QLLMS diff_Rot_QLLMS diff_T_QLRoot  diff_Rot_QLRoot diff_Rot_RA4  diff_T_RA4;
% figure(30)
% plot(mean(diff_Rot_RA4,2),'x')
% hold 
% plot(mean(diff_Rot_QLLMS,2),'*')
% hold on
% plot(mean(diff_Rot_QLRoot,2),'d')
% xlabel('Noise level')
% ylabel('Median Rotation error')
% axis([-8,10,0,1])
% legend('RA','LMS','Root')
% figure(32)
% plot(mean(diff_T_RA4,2),'x')
% hold 
% plot(mean(diff_T_QLLMS,2),'*')
% hold on
% plot(mean(diff_T_QLRoot,2),'d')
% xlabel('Noise level')
% ylabel('Median Translation error')
% legend('RA','LMS','Root')
% axis([0,10,0,50])
% 
% load data_PNP_C1
% figure(30)
% plot(median(diff_Rot_epnp,2),'o')
% hold 
% plot(median(diff_Rot_Paul,2),'*')
% xlabel('Noise level')
% ylabel('Median Rotation error')
% % legend('Sigle','LMS','Root')
% figure(32)
% plot(median(diff_T_epnp,2),'o')
% hold 
% plot(median(diff_T_Paul,2),'*')
% xlabel('Noise level')
% ylabel('Median Translation error')
% % legend('Sigle','LMS','Root')


load data_PNP_C2_100
figure(30)
plot(median(diff_Rot_epnp,2),'d')
hold on
plot(median(diff_Rot_Paul,2),'*')
hold on
plot(median(diff_Rot_mean50,2),'x')
hold on
plot(median(diff_Rot_mean100,2),'.')
axis([0,10,0,1])
xlabel('Noise level')
ylabel('Median Rotation error')
% legend('Lepetit','Fiore','RA50','RA100')
% legend('Lepetit','Fiore','RA')
figure(32)
plot(median(diff_T_epnp,2),'d')
hold on
plot(median(diff_T_Paul,2),'*')
hold on
plot(median(diff_T_mean50,2),'x')
hold on
plot(median(diff_T_mean100,2),'.')
axis([0,10,0,5])
xlabel('Noise level')
ylabel('Median Translation error')
% legend('Lepetit','Fiore','RA50','RA100')
% legend('Lepetit','Fiore','RA')

% load data_PNP_C2_100
figure(31)
plot(mean(diff_Rot_epnp,2),'d')
hold on
plot(mean(diff_Rot_Paul,2),'*')
hold on
plot(mean(diff_Rot_mean50,2),'x')
hold on
plot(mean(diff_Rot_mean100,2),'.')
axis([0,10,0,1])
xlabel('Noise level')
ylabel('Mean Rotation error')
% legend('Lepetit','Fiore','RA50','RA100')
% legend('Lepetit','Fiore','RA')
figure(33)
plot(mean(diff_T_epnp,2),'d')
hold on
plot(mean(diff_T_Paul,2),'*')
hold on
plot(mean(diff_T_mean50,2),'x')
hold on
plot(mean(diff_T_mean100,2),'.')
axis([0,10,0,5])
xlabel('Noise level')
ylabel('Mean Translation error')
% legend('Lepetit','Fiore','RA50','RA100')
pause
% legend('Lepetit','Fiore','RA')
% 
