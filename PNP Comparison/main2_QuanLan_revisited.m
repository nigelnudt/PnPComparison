clear all; close all;
%main1 主要为了比较Quan&Lan方法在4点以及5点时究竟是否比解方程容易
addpath 'C:\Users\nigel\Documents\Multiprecision Computing Toolbox'
addpath Laurent_Kneip_CVPR_2011\p3p_code_final;
addpath Lu_Hager_Mjolsness_TPAMI_2000;
addpath Vincent_lepetit_EPnP_IJCV_2008
addpath Vincent_lepetit_EPnP_IJCV_2008\data;
addpath Vincent_lepetit_EPnP_IJCV_2008\error;
addpath Vincent_lepetit_EPnP_IJCV_2008\EPnP;
% addpath D:\计算机视觉测量\视觉测量（程序）\PNP Comparison\Laurent_Kneip_CVPR_2011\p3p_code_final
num_points=50;
f_camera=800;
% std_noise=10/f_camera;%得到的误差属于
%1.-Generate simulated input data------------------------------------------

diff_Rot_Kneip=zeros(10,100);
% diff_Rot_Fischler=zeros(10,100);
diff_Rot_QLRoot=zeros(10,100);
diff_Rot_QLLMS=zeros(10,100);
diff_Rot_epnp=zeros(10,100);
diff_Rot_RA4=zeros(10,100);

%         diff_T_Fischler(noise_level,iterative_index)=norm(T_fishler-T_real);        
diff_T_Kneip=zeros(10,100);
% diff_T_Fischler=zeros(10,100);
diff_T_QLRoot=zeros(10,100);
diff_T_QLLMS=zeros(10,100);
diff_T_epnp=zeros(10,100);
diff_T_RA4=zeros(10,100);
% for noise_level=1:10
for noise_level=10:-1:1
%     std_noise=noise_level*0;
    T_LongQun_root_sum=0;
    R_LongQun_root_sum=0;
    for iterative_index=1:100
        load_points=0;
        if ~load_points
            n=num_points; %number of points
%             std_noise=1; %noise in the measurements (in pixels)
            [A,point,Rt,T_real]=generate_noisy_input_data(n,noise_level);
            save('data\input_data_noise.mat','A','point','Rt');
        else
            load('data\input_data_noise.mat','A','point','Rt');
            n=size(point,2);
            draw_noisy_input_data(point);
        end
        
        %2.-Inputs format--------------------------------
        x3d=zeros(n,4);
        x2d=zeros(n,3);
        A=A(:,1:3);
        for i=1:n
            x3d_h(i,:)=[point(i).Xworld',1];
            x2d_h(i,:)=[point(i).Ximg(1:2)',1];
            
            %world and camera coordinates
            X3d_world(i,:)=point(i).Xworld';
            X3d_cam(i,:)=point(i).Xcam';
            X_distance_cam(i)=norm(point(i).Xcam');
        end
        
        Xw=x3d_h(:,1:3);
        U=x2d_h(:,1:2);
        x2d_h_normlized=(A\x2d_h')';
        Rot_origin=Rt(1:3,1:3);
        [theta_origin, vector_origin]=rodrigues_rot2vetor(Rot_origin) ;
        p_origin=theta_origin* vector_origin;
        
        [Repnp,Tepep,Xc,sol]=efficient_pnp(x3d_h(1:4,:),x2d_h(1:4,:),A);%Lepetit,2006,EPnP
        f1=x2d_h_normlized(1,:)/norm(x2d_h_normlized(1,:));
        f2=x2d_h_normlized(2,:)/norm(x2d_h_normlized(2,:));
        f3=x2d_h_normlized(3,:)/norm(x2d_h_normlized(3,:));
        f4=x2d_h_normlized(4,:)/norm(x2d_h_normlized(4,:));
        poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
        poses2 = p3p( X3d_world([1,2,4],1:3)', [f1',f2',f4'] );
        poses3 = p3p( X3d_world([1,3,4],1:3)', [f1',f3',f4'] );
        poses4 = p3p( X3d_world([2,3,4],1:3)', [f2',f3',f4'] );
        
        [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
        [R_Kneip2,T_Kneip2]=solution_select_kneip(poses2,X3d_world(3,1:3),x2d_h_normlized(3,:));
        [R_Kneip3,T_Kneip3]=solution_select_kneip(poses3,X3d_world(2,1:3),x2d_h_normlized(2,:));
        [R_Kneip4,T_Kneip4]=solution_select_kneip(poses4,X3d_world(1,1:3),x2d_h_normlized(1,:));
        [R_LongQun_root,T_LongQun_root,Xc3,R_LongQun_LMS,T_LongQun_LMS,Xc4]=LongQuan3Point4PointPnP(x3d_h,x2d_h_normlized);
        
        
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         tic
%         for index=1:100
        
        R_LongQun_root_sum=R_LongQun_root_sum+R_LongQun_root;
        T_LongQun_root_sum=T_LongQun_root_sum+T_LongQun_root;
        
        
        [theta_epnp, vector_epnp]=rodrigues_rot2vetor(Repnp) ;            
        p_epnp=theta_epnp*vector_epnp;
        diff_Rot_epnp(noise_level,iterative_index)=min(norm(p_origin-p_epnp),norm(-p_origin-p_epnp));
        
        
        [theta_Kneip, vector_Kneip]=rodrigues_rot2vetor(R_Kneip) ;            
        p_Kneip=theta_Kneip*vector_Kneip;
        diff_Rot_Kneip(noise_level,iterative_index)=min(norm(p_origin-p_Kneip),norm(-p_origin-p_Kneip));
        
        [theta_Kneip, vector_Kneip]=rodrigues_rot2vetor(R_Kneip) ;            
        p_Kneip=theta_Kneip*vector_Kneip;
        [theta_Kneip2, vector_Kneip2]=rodrigues_rot2vetor(R_Kneip2) ;            
        p_Kneip2=theta_Kneip2*vector_Kneip2;
        [theta_Kneip3, vector_Kneip3]=rodrigues_rot2vetor(R_Kneip3) ;            
        p_Kneip3=theta_Kneip*vector_Kneip;
        [theta_Kneip4, vector_Kneip4]=rodrigues_rot2vetor(R_Kneip4) ;            
        p_Kneip4=theta_Kneip4*vector_Kneip4;
        R_rot_average=(p_Kneip4+p_Kneip3+p_Kneip2+p_Kneip)/4;
        T_rot_average=(T_Kneip+T_Kneip2+T_Kneip3+T_Kneip4)/4;
%         [R_rot_average,T_rot_average]=rotation_average_kneip(poses,poses2,poses3,poses4);
        
        
        [theta__QLRoot, vector_QLRoot]=rodrigues_rot2vetor(R_LongQun_root) ;            
        p_QLRoot=theta__QLRoot*vector_QLRoot;
        diff_Rot_QLRoot(noise_level,iterative_index)=min(norm(p_origin-p_QLRoot),norm(-p_origin-p_QLRoot));        
        
        [theta_QLLMS, vector_QLLMS]=rodrigues_rot2vetor(R_LongQun_LMS) ;            
        p_QLLMS=theta_QLLMS*vector_QLLMS;
        diff_Rot_QLLMS(noise_level,iterative_index)=min(norm(p_origin-p_QLLMS),norm(-p_origin-p_QLLMS));
        
%         [theta_RA4, vector_RA4]=rodrigues_rot2vetor(R_rot_average) ;            
%         p_RA4=theta_RA4*vector_RA4;
        diff_Rot_RA4(noise_level,iterative_index)=min(norm(p_origin-R_rot_average),norm(-p_origin-R_rot_average));
        
%         diff_Rot_epnp(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Repnp-eye(3))));
% %         diff_Rot_Paul(noise_level,iterative_index)=sum(sum(abs(Rot_origin*RPaul-eye(3))));
%         diff_Rot_kneip(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_Kneip-eye(3))));
% %         diff_Rot_Fischler(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_fischler-eye(3))));
%         diff_Rot_QuanLan_root(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_LongQun_root-eye(3))));
%         diff_Rot_QuanLan_LMS(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_LongQun_LMS-eye(3))));
%         diff_Rot_4mean(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_rot_average-eye(3))));
        
        diff_T_epnp(noise_level,iterative_index)=norm(Tepep-T_real);   
        diff_T_Kneip(noise_level,iterative_index)=norm(T_Kneip-T_real);
%         diff_T_Fischler(noise_level,iterative_index)=norm(T_fishler-T_real);        
        diff_T_QLRoot(noise_level,iterative_index)=norm(T_LongQun_root-T_real); 
        diff_T_QLLMS(noise_level,iterative_index)=norm(T_LongQun_LMS-T_real);
        diff_T_RA4(noise_level,iterative_index)=norm(T_rot_average-T_real);
        
%         diff_T_3p(noise_level,iterative_index)=norm(tp3-Rp3'*T_real');
%         diff_T_4p(noise_level,iterative_index)=norm(tp4-Rp4'*T_real');
%         diff_T_epnp(noise_level,iterative_index)=norm(Tp-T_real');
%         diff_T_Paul(noise_level,iterative_index)=norm(TPaul-RPaul'*T_real');
        disp(num2str(iterative_index));
    end
end
figure(2)
plot(diff_Rot_Kneip(3,:),'.');
axis([0,100,0,2])
title('Kneip'); 
ylabel('Rotation error')
% hold on%         diff_Rot_3p(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Rp3-eye(3))));
%         diff_Rot_4p(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Rp4-eye(3))));
% figure(3)
% plot(diff_Rot_Fischler(10,:),'.');
% axis([0,100,0,5])
% title('Fischler & Bolles'); 
% ylabel('Rotation error')
figure(4)
% hold on
plot(diff_Rot_QLRoot(3,:),'.');
a11=sum(diff_Rot_QLRoot(3,:));
axis([0,100,0,2])
title('Quan & Lan root'); 
ylabel('Rotation error')
% hold on
figure(5)
plot(diff_Rot_QLLMS(3,:),'.');
axis([0,100,0,2])
title('Quan & Lan LMS'); 
ylabel('Rotation error')
figure(6)
plot(diff_Rot_epnp(3,:),'.');
axis([0,100,0,2])
title('Lepetit'); 
ylabel('Rotation error')
figure(7)
plot(diff_Rot_RA4(3,:),'.');
axis([0,100,0,2])
a12=sum(diff_Rot_RA4(3,:));
title('Rotation Averaging'); 
ylabel('Rotation error')
% hold on
figure(12)
plot(diff_T_Kneip(3,:),'.');
axis([0,100,0,20])
title('Kneip'); 
ylabel('Translation error')
% % hold on
% figure(13)
% plot(diff_T_Fischler(10,:),'.');
% axis([0,100,0,20])
% title('Fischler & Bolles'); 
% ylabel('Translation error')
figure(14)
% hold on
plot(diff_T_QLRoot(3,:),'.');
axis([0,100,0,50])
title('Quan & Lan root'); 
ylabel('Translation error')
% hold on
figure(15)
plot(diff_T_QLLMS(3,:),'.');
axis([0,100,0,50])
title('Quan & Lan LMS'); 
ylabel('Translation error')
pause(0.1)
figure(16)
plot(diff_T_epnp(3,:),'.');
axis([0,100,0,50])
title('Lepetit'); 
ylabel('Translation error')
figure(17)
plot(diff_T_RA4(3,:),'.');
axis([0,100,0,50])
title('Rotation Averaging'); 
ylabel('Translation error')
pause(0.1)
figure(30)
plot(mean(diff_T_QLLMS,2),'*')
hold on
plot(mean(diff_T_QLRoot,2),'.')
% hold on
% plot(mean(diff_T_epnp,2),'.')
hold on
plot(mean(diff_T_RA4,2),'x')
% hold on
% plot(mean(diff_T_Kneip,2),'o')
ylabel('Average Translation error')
legend('LMS','Root','RA4')
axis([1,10,0,50])
% legend('LMS','Root','Lepetit','RA4','Kneip')
figure(32)
plot(mean(diff_Rot_QLLMS,2),'*')
hold on
plot(mean(diff_Rot_QLRoot,2),'.')
% hold on
% plot(mean(diff_Rot_epnp,2),'.')
hold on
plot(mean(diff_Rot_RA4,2),'x')
axis([1,10,0,1])
% hold on
% plot(mean(diff_Rot_Kneip,2),'o')
ylabel('Average Rotation error')
legend('LMS','Root','RA4')
save data_QuanLan_C3 diff_T_epnp diff_Rot_epnp diff_T_QLLMS diff_Rot_QLLMS diff_T_QLRoot  diff_Rot_QLRoot diff_Rot_RA4  diff_T_RA4;













