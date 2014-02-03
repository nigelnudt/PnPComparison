clear all; close all;
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
std_noise=10/f_camera;%得到的误差属于
%1.-Generate simulated input data------------------------------------------
for iterative_index=1:100
    load_points=0;
    if ~load_points
        n=num_points; %number of points
        std_noise=1; %noise in the measurements (in pixels)
%         [A,point,Rt]=generate_noisy_input_data(n,std_noise);
        [A,point,Rt,centroid]=generate_noisy_input_data(n,std_noise)
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
    
    [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
    [Rp,Tp,Xc,sol]=efficient_pnp(x3d_h(1:50,:),x2d_h(1:50,:),A);%Lepetit,2006,EPnP
    [RPaul,TPaul]=Paul_efficient_two_step_PNP(x3d_h(1:50,:),x2d_h_normlized(1:50,:));
    f1=x2d_h_normlized(1,:)/norm(x2d_h_normlized(1,:));
    f2=x2d_h_normlized(2,:)/norm(x2d_h_normlized(2,:));
    f3=x2d_h_normlized(3,:)/norm(x2d_h_normlized(3,:));
    poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
    poses2 = p3p( X3d_world([2,3,1],1:3)', [f2',f3',f1'] );
    poses3 = p3p( X3d_world([3,1,2],1:3)', [f3',f1',f2'] );
%     [Rp3,tp3,Xc3,Rp4,tp4,Xc4]=LongQuan3Point4PointPnP(x3d_h,x2d_h_normlized);
    [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3]=LongQuan3Point(x3d_h,x2d_h_normlized);
    
    
    
    %     diff_kneip1=sum(sum(abs(Rot_origin*poses(:,2:4)'-eye(3))));
    %     diff_kneip2=sum(sum(abs(Rot_origin*poses(:,6:8)'-eye(3))));
    %     diff_kneip3=sum(sum(abs(Rot_origin*poses(:,10:12)'-eye(3))));
    %     diff_kneip4=sum(sum(abs(Rot_origin*poses(:,14:16)'-eye(3))));
    %     [diff_Rot_kneip(iterative_index),kneip_index]=min([diff_kneip1,diff_kneip2,diff_kneip3,diff_kneip4]);
    %     diff_Rot_3p(iterative_index)=sum(sum(abs(Rot_origin*Rp3-eye(3))));
    %     diff_Rot_4p(iterative_index)=sum(sum(abs(Rot_origin*Rp4-eye(3))));
    %     diff_Rot_epnp(iterative_index)=sum(sum(abs(Rot_origin*Rp-eye(3))));
    %     diff_Rot_Paul(iterative_index)=sum(sum(abs(Rot_origin*RPaul-eye(3))));
    %     diff_Rot_Fischler(iterative_index)=sum(sum(abs(Rot_origin*RPaul-eye(3))));
    [theta_kneip1, vector_kneip1]=rodrigues_rot2vetor(poses(:,2:4)) ;
    p_kneip1=theta_kneip1*vector_kneip1;
    diff_kneip1=min(norm(p_origin-p_kneip1),norm(-p_origin-p_kneip1));
    
    
    [theta_kneip2, vector_kneip2]=rodrigues_rot2vetor(poses(:,6:8));
    p_kneip2=theta_kneip2*vector_kneip2;
    diff_kneip2=min(norm(p_origin-p_kneip2),norm(-p_origin-p_kneip2));
    
    [theta_kneip3, vector_kneip3]=rodrigues_rot2vetor(poses(:,10:12)) ;
    p_kneip3=theta_kneip3*vector_kneip3;
    diff_kneip3=min(norm(p_origin-p_kneip3),norm(-p_origin-p_kneip3));
    
    [theta_kneip4, vector_kneip4]=rodrigues_rot2vetor(poses(:,14:16)) ;
    p_kneip4=theta_kneip4*vector_kneip4;
    diff_kneip4=min(norm(p_origin-p_kneip4),norm(-p_origin-p_kneip4));
%     for kneip_index=1:4
%         p4_camera=Rp3*((x3d_h(4,1:3))'+tp3);
%         p4_error=norm([p4_camera(1)/p4_camera(3),p4_camera(2)/p4_camera(3)]-x2d_h_normlized(4,1:2));
%         if p4_error<verify_error
%             R_LongQun=Rp3;
%             T_LongQun=Rp3*tp3;
%             verify_error=p4_error;
%         end
%     end
%     residual_total(1,solution_index)=residual2+residual3;
    [diff_Rot_kneip(iterative_index),kneip_index]=min([diff_kneip1,diff_kneip2,diff_kneip3,diff_kneip4]);
    
    [theta_3p, vector_3p]=rodrigues_rot2vetor(R_LongQun) ;
    p_3p=theta_3p*vector_3p;
    diff_Rot_3p(iterative_index)=min(norm(p_origin-p_3p),norm(-p_origin-p_3p));
    
    [theta_4p, vector_4p]=rodrigues_rot2vetor(R_LongQun) ;
    p_4p=theta_4p*vector_4p;
    diff_Rot_4p(iterative_index)=min(norm(p_origin-p_4p),norm(-p_origin-p_4p));
    
    [theta_epnp, vector_epnp]=rodrigues_rot2vetor(Rp) ;
    p_epnp=theta_epnp*vector_epnp;
    diff_Rot_epnp(iterative_index)=min(norm(p_origin-p_epnp),norm(-p_origin-p_epnp));
    
    [theta_Paul, vector_Paul]=rodrigues_rot2vetor(RPaul) ;
    p_Paul=theta_Paul*vector_Paul;
    diff_Rot_Paul(iterative_index)=min(norm(p_origin-p_Paul),norm(-p_origin-p_Paul));
    
    [theta_Fischler, vector_Fischler]=rodrigues_rot2vetor(R_fischler) ;
    p_Fischler=theta_Fischler*vector_Fischler;
    diff_Rot_Fischler(iterative_index)=min(norm(p_origin-p_Fischler),norm(-p_origin-p_Fischler));
    
%     diff_Rot_3p(iterative_index)=sum(sum(abs(Rot_origin*Rp3-eye(3))));
%     diff_Rot_4p(iterative_index)=sum(sum(abs(Rot_origin*Rp4-eye(3))));
%     diff_Rot_epnp(iterative_index)=sum(sum(abs(Rot_origin*Rp-eye(3))));
%     diff_Rot_Paul(iterative_index)=sum(sum(abs(Rot_origin*RPaul-eye(3))));
%     diff_Rot_Fischler(iterative_index)=sum(sum(abs(Rot_origin*R_fischler-eye(3))));
    diff_T_kneip(iterative_index)=norm(poses(:,4*kneip_index-3)+poses(:,(4*kneip_index-2):(4*kneip_index))*centroid);
    diff_T_3p(iterative_index)=norm(T_LongQun-centroid);
    diff_T_4p(iterative_index)=norm(T_LongQun-centroid);
    diff_T_epnp(iterative_index)=norm(Tp-centroid);
    diff_T_Paul(iterative_index)=norm(TPaul-centroid);
    diff_T_Fischler(iterative_index)=norm(T_fishler-centroid);
    disp(num2str(iterative_index));
end
figure(2)
plot(diff_Rot_kneip,'*');
axis([0,100,0,0.5])
% hold on
figure(3)
plot(diff_Rot_3p,'*');
axis([0,100,0,0.5])
figure(4)
% hold on
plot(diff_Rot_4p,'*');
axis([0,100,0,0.5])
% hold on
figure(5)
plot(diff_Rot_epnp,'*');
axis([0,100,0,0.5])
figure(6)
plot(diff_Rot_Paul,'*');
axis([0,100,0,0.5])
figure(7)
plot(diff_Rot_Fischler,'*');
axis([0,100,0,0.5])
% hold on
figure(12)
plot(diff_T_kneip,'*');
axis([0,100,0,1])
% hold on
figure(13)
plot(diff_T_3p,'*');
axis([0,100,0,1])
figure(14)
% hold on
plot(diff_T_4p,'*');
axis([0,100,0,1])
% hold on
figure(15)
plot(diff_T_epnp,'*');
axis([0,100,0,1])
pause(0.1)
figure(16)
plot(diff_T_Paul,'*');
axis([0,100,0,1])
figure(17)
plot(diff_T_Fischler,'*');
axis([0,100,0,1])
pause(0.1)
save data_general_depth_SFOV diff_Rot_kneip diff_Rot_3p diff_Rot_4p diff_Rot_epnp diff_Rot_Paul diff_Rot_Fischler diff_T_kneip diff_T_3p diff_T_4p diff_T_epnp diff_T_Paul diff_T_Fischler













