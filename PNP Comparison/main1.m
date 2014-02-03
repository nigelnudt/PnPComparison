clear all; close all;
%main1 主要为了比较P3P，包括Fishler&bolles，Quan&Lan，Kneip
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

diff_Rot_kneip=zeros(10,100);
diff_Rot_Fischler=zeros(10,100);
diff_Rot_QuanLan=zeros(10,100);
diff_Rot_QuanLan_mp=zeros(10,100);
diff_Rot_epnp=zeros(10,100);

diff_T_kneip=zeros(10,100);
diff_T_Fischler=zeros(10,100);
diff_T_QuanLan=zeros(10,100);
diff_T_QuanLan_mp=zeros(10,100);
diff_T_epnp=zeros(10,100);

% for noise_level=1:10
for noise_level=10:-1:1
    std_noise=noise_level;
    for iterative_index=1:100
        load_points=0;
        if ~load_points
            n=num_points; %number of points
%             std_noise=1; %noise in the measurements (in pixels)
            [A,point,Rt,T_real]=generate_noisy_input_data(n,std_noise*0.3);
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
        p_origin=theta_origin*vector_origin;
        
        [Repnp,Tepep,Xc,sol]=efficient_pnp(x3d_h(1:50,:),x2d_h(1:50,:),A);%Lepetit,2006,EPnP
%         [RPaul,TPaul]=Paul_efficient_two_step_PNP(x3d_h(1:50,:),x2d_h_normlized(1:50,:));
        f1=x2d_h_normlized(1,:)/norm(x2d_h_normlized(1,:));
        f2=x2d_h_normlized(2,:)/norm(x2d_h_normlized(2,:));
        f3=x2d_h_normlized(3,:)/norm(x2d_h_normlized(3,:));
%         tic
%         for i=1:1000
        poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
        [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
%         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
%         
%         end
%         toc
%         tic
%         for i=1:1000
        [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
%         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
        
%         end
%         toc
%         poses2 = p3p( X3d_world([2,3,1],1:3)', [f2',f3',f1'] );
%         poses3 = p3p( X3d_world([3,1,2],1:3)', [f3',f1',f2'] );
%         tic
%         for i=1:1000
        [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         end
%         toc
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
%         [Rp3,tp3,Xc3,Rp4,tp4,Xc4]=LongQuan3Point4PointPnP(x3d_h,x2d_h_normlized);
        
    [theta_Kneip, vector_Kneip]=rodrigues_rot2vetor(R_Kneip) ;
    p_Kneip=theta_Kneip*vector_Kneip;
    diff_Rot_kneip(noise_level,iterative_index)=min(norm(p_origin-p_Kneip),norm(-p_origin-p_Kneip));
 
    
    [theta_QL, vector_QL]=rodrigues_rot2vetor(R_LongQun) ;
    p_QL=theta_QL*vector_QL;
    diff_Rot_QuanLan(noise_level,iterative_index)=min(norm(p_origin-p_QL),norm(-p_origin-p_QL));
    
    [theta_QL_mp, vector_QL_mp]=rodrigues_rot2vetor(R_LongQun_mp) ;
    p_QL_mp=theta_QL_mp*vector_QL_mp;
    diff_Rot_QuanLan_mp(noise_level,iterative_index)=min(norm(p_origin-p_QL_mp),norm(-p_origin-p_QL_mp));
    
%     [theta_epnp, vector_epnp]=rodrigues_rot2vetor(Rp) ;
%     p_epnp=theta_epnp*vector_epnp;
%     diff_Rot_epnp(iterative_index)=min(norm(p_origin-p_epnp),norm(-p_origin-p_epnp));
    
%     [theta_Paul, vector_Paul]=rodrigues_rot2vetor(RPaul) ;
%     p_Paul=theta_Paul*vector_Paul;
%     diff_Rot_Paul(iterative_index)=min(norm(p_origin-p_Paul),norm(-p_origin-p_Paul));
    
    [theta_Fischler, vector_Fischler]=rodrigues_rot2vetor(R_fischler) ;
    p_Fischler=theta_Fischler*vector_Fischler;
    diff_Rot_Fischler(noise_level,iterative_index)=min(norm(p_origin-p_Fischler),norm(-p_origin-p_Fischler));
% 
%         diff_Rot_epnp(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Repnp-eye(3))));
% %         diff_Rot_Paul(noise_level,iterative_index)=sum(sum(abs(Rot_origin*RPaul-eye(3))));
%         diff_Rot_kneip(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_Kneip-eye(3))));
%         diff_Rot_Fischler(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_fischler-eye(3))));
%         diff_Rot_QuanLan(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_LongQun-eye(3))));
%         diff_Rot_QuanLan_mp(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_LongQun_mp-eye(3))));
        
        diff_T_kneip(noise_level,iterative_index)=norm(T_Kneip-T_real);
        diff_T_Fischler(noise_level,iterative_index)=norm(T_fishler-T_real);        
        diff_T_QuanLan(noise_level,iterative_index)=norm(T_LongQun-T_real);
        diff_T_QuanLan_mp(noise_level,iterative_index)=norm(T_LongQun_mp-T_real);
        
        
%         diff_T_3p(noise_level,iterative_index)=norm(tp3-Rp3'*T_real');
%         diff_T_4p(noise_level,iterative_index)=norm(tp4-Rp4'*T_real');
%         diff_T_epnp(noise_level,iterative_index)=norm(Tp-T_real');
%         diff_T_Paul(noise_level,iterative_index)=norm(TPaul-RPaul'*T_real');
        disp(num2str(iterative_index));
    end
end
figure(2)
plot(diff_Rot_kneip(10,:),'.');
axis([0,100,0,2])
title('Kneip'); 
ylabel('Rotation error')
% hold on%         diff_Rot_3p(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Rp3-eye(3))));
%         diff_Rot_4p(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Rp4-eye(3))));
figure(3)
plot(diff_Rot_Fischler(10,:),'.');
axis([0,100,0,2])
title('Fischler & Bolles'); 
ylabel('Rotation error')
figure(4)
% hold on
plot(diff_Rot_QuanLan(10,:),'.');
axis([0,100,0,2])
title('Quan & Lan'); 
ylabel('Rotation error')
% hold on
figure(5)
plot(diff_Rot_QuanLan_mp(10,:),'.');
axis([0,100,0,2])
title('Quan & Lan multi-precision'); 
ylabel('Rotation error')
figure(6)
plot(diff_Rot_epnp(10,:),'.');
axis([0,100,0,2])
title('Lepetit'); 
ylabel('Rotation error')
% hold on
figure(12)
plot(diff_T_kneip(10,:),'.');
axis([0,100,0,10])
title('Kneip'); 
ylabel('Translation error')
% hold on
figure(13)
plot(diff_T_Fischler(10,:),'.');
axis([0,100,0,10])
title('Fischler & Bolles'); 
ylabel('Translation error')
figure(14)
% hold on
plot(diff_T_QuanLan(10,:),'.');
axis([0,100,0,10])
title('Quan & Lan'); 
ylabel('Translation error')
% hold on
figure(15)
plot(diff_T_QuanLan_mp(10,:),'.');
axis([0,100,0,10])
title('Quan & Lan multi-precision'); 
ylabel('Translation error')
pause(0.1)
figure(16)
plot(diff_T_epnp(10,:),'.');
axis([0,100,0,10])
title('Lepetit'); 
ylabel('Translation error')
pause(0.1)
save data_P3P_C1 diff_Rot_kneip diff_Rot_Fischler diff_Rot_QuanLan diff_Rot_QuanLan_mp diff_T_QuanLan diff_T_QuanLan_mp diff_T_Fischler diff_T_kneip;













