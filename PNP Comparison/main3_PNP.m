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
diff_Rot_Paul=zeros(10,100);
diff_Rot_mean50=zeros(10,100);
diff_Rot_mean100=zeros(10,100);

diff_T_kneip=zeros(10,100);
diff_T_Fischler=zeros(10,100);
diff_T_QuanLan=zeros(10,100);
diff_T_QuanLan_mp=zeros(10,100);
diff_T_epnp=zeros(10,100);
diff_T_Paul=zeros(10,100);
diff_T_mean50=zeros(10,100);
diff_T_mean100=zeros(10,100);
% for noise_level=1:10
isnot_initial_mean=zeros(1,10);
for noise_level=10:-1:1
    std_noise=noise_level;
    for iterative_index=1:100
        load_points=0;
        if ~load_points
            n=num_points; %number of points
            %             std_noise=1; %noise in the measurements (in pixels)
            [A,point,Rt,T_real]=generate_noisy_input_data(n,std_noise*1);
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
        for pindex=1:n
            x2d_h_vector(pindex,:)=x2d_h_normlized(pindex,:)/norm(x2d_h_normlized(pindex,:));
        end
        Rot_origin=Rt(1:3,1:3);
        [theta11_origin, vector_y_origin]=rodrigues_rot2vetor(Rot_origin) ;
        p_origin=theta11_origin* vector_y_origin;
        %
        %         Rot_recover=eye(3)+sin(theta)*v_asym_matrix+(1-cos(theta))*v_asym_matrix*v_asym_matrix;
        [Repnp,Tepep,Xc,sol]=efficient_pnp(x3d_h(1:50,:),x2d_h(1:50,:),A);%Lepetit,2006,EPnP
        [RPaul,TPaul]=Paul_efficient_two_step_PNP(x3d_h(1:50,:),x2d_h_normlized(1:50,:));
        
        is_find_proper_set=0;
        for index=1:46
            if is_find_proper_set==0
                a=abs(det(X3d_world(1:3,1:3)));
                b=abs(det(X3d_world([1,2,4],1:3)));
                c=abs(det(X3d_world([1,3,4],1:3)));
                d=abs(det(X3d_world([2,3,4],1:3)));
                if a>1 && b>1 && c>1 && d>1
                    poses = p3p( X3d_world(1:3,1:3)', [x2d_h_vector(1,:)',x2d_h_vector(2,:)',x2d_h_vector(3,:)'] );
                    poses2 = p3p( X3d_world([1,2,4],1:3)', [x2d_h_vector(1,:)',x2d_h_vector(2,:)',x2d_h_vector(4,:)'] );
                    poses3 = p3p( X3d_world([1,3,4],1:3)', [x2d_h_vector(1,:)',x2d_h_vector(3,:)',x2d_h_vector(4,:)'] );
                    poses4 = p3p( X3d_world([2,3,4],1:3)', [x2d_h_vector(2,:)',x2d_h_vector(3,:)',x2d_h_vector(4,:)'] );
                    [R_rot_average,T_rot_average]=rotation_average_kneip(poses,poses2,poses3,poses4);
                    is_find_proper_set=1;
                end
            end
        end
        if is_find_proper_set==0
            R_rot_average=Repnp;
            T_rot_average=Tepep;
            isnot_initial_mean(noise_level)=isnot_initial_mean(noise_level)+1;
        end
        %         tic
        for index=1:50
            all_set=[1:50];
            index_p1=ceil(50*rand);
            all_set=setdiff(all_set,[index_p1]);
            index_p2=all_set(ceil(49*rand));
            all_set=setdiff(all_set,[index_p2]);
            index_p3=all_set(ceil(48*rand));
            coordinate_this_set=[ X3d_world([index_p1,index_p2,index_p3],1:3)];
            if abs(det(coordinate_this_set))>2%只有不在同一条先线上才进行计算。
                pose5 = p3p( X3d_world([index_p1,index_p2,index_p3],1:3)', [x2d_h_vector(index_p1,:)',x2d_h_vector(index_p2,:)',x2d_h_vector(index_p3,:)']);
                weight=1/(index+1);
                [R_rot_average,T_rot_average]=rotation_average_weighted(R_rot_average,T_rot_average,pose5,weight);
            end
%             diff_Rot_Mean(index)=sum(sum(abs(Rot_origin*R_rot_average'-eye(3))));
%             diff_T_Mean(index)=norm(T_rot_average-T_real);
        end
        R_rot_average50=R_rot_average;
        T_rot_average50=T_rot_average;
        for index=1:50
            all_set=[1:50];
            index_p1=ceil(50*rand);
            all_set=setdiff(all_set,[index_p1]);
            index_p2=all_set(ceil(49*rand));
            all_set=setdiff(all_set,[index_p2]);
            index_p3=all_set(ceil(48*rand));
            coordinate_this_set=[ X3d_world([index_p1,index_p2,index_p3],1:3)];
            if abs(det(coordinate_this_set))>2%只有不在同一条先线上才进行计算。
                pose5 = p3p( X3d_world([index_p1,index_p2,index_p3],1:3)', [x2d_h_vector(index_p1,:)',x2d_h_vector(index_p2,:)',x2d_h_vector(index_p3,:)']);
                weight=1/(index+1);
                [R_rot_average,T_rot_average]=rotation_average_weighted(R_rot_average,T_rot_average,pose5,weight);
            end
%             diff_Rot_Mean(index)=sum(sum(abs(Rot_origin*R_rot_average'-eye(3))));
%             diff_T_Mean(index)=norm(T_rot_average-T_real);
        end
        %         toc
%         pause(0.01)
        %         tic
        %         [RPaul,TPaul]=Paul_efficient_two_step_PNP(x3d_h(1:50,:),x2d_h_normlized(1:50,:));
        %         f1=x2d_h_normlized(1,:)/norm(x2d_h_normlized(1,:));
        %         f2=x2d_h_normlized(2,:)/norm(x2d_h_normlized(2,:));
        %         f3=x2d_h_normlized(3,:)/norm(x2d_h_normlized(3,:));
        % %         tic
        % %         for i=1:1000
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
        %         poses = p3p( X3d_world(1:3,1:3)', [f1',f2',f3'] );
        %         [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world(4,1:3),x2d_h_normlized(4,:));
        %
        %         end
        %         toc
        %         tic
        %         for i=1:1000
        %         [R_fischler,T_fishler]=P3P_fischler_bolles(x3d_h,x2d_h_normlized);
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
        %         [R_LongQun,T_LongQun,R_LongQun_mp,T_LongQun_mp,Xc3_LongQuan]=LongQuan3Point(x3d_h,x2d_h_normlized);
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
%         p_origin=theta11_origin* vector_y_origin;
        [theta_epnp, vector_epnp]=rodrigues_rot2vetor(Repnp) ;            
        p_epnp=theta_epnp*vector_epnp;
        diff_Rot_epnp(noise_level,iterative_index)=min(norm(p_origin-p_epnp),norm(-p_origin-p_epnp));
        [theta_Fiore, vector_Fiore]=rodrigues_rot2vetor(RPaul) ;            
        p_Fiore=theta_Fiore*vector_Fiore;
        diff_Rot_Paul(noise_level,iterative_index)=min(norm(p_origin-p_Fiore),norm(-p_origin-p_Fiore));
        [theta_mean50, vector_mean50]=rodrigues_rot2vetor(R_rot_average50) ;            
        p_mean50=theta_mean50*vector_mean50;
        diff_Rot_mean50(noise_level,iterative_index)=min(norm(p_origin-p_mean50),norm(-p_origin-p_mean50));
        [theta_mean, vector_mean]=rodrigues_rot2vetor(R_rot_average) ;            
        p_mean=theta_mean*vector_mean;
        diff_Rot_mean100(noise_level,iterative_index)=min(norm(p_origin-p_mean),norm(-p_origin-p_mean));
        
        %         diff_Rot_kneip(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_Kneip-eye(3))));
        %         diff_Rot_Fischler(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_fischler-eye(3))));
        %         diff_Rot_QuanLan(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_LongQun-eye(3))));
        %         diff_Rot_QuanLan_mp(noise_level,iterative_index)=sum(sum(abs(Rot_origin*R_LongQun_mp-eye(3))));
        
        diff_T_epnp(noise_level,iterative_index)=norm(Tepep-T_real);
        diff_T_Paul(noise_level,iterative_index)=norm(TPaul-T_real);
        diff_T_mean50(noise_level,iterative_index)=norm(T_rot_average50-T_real);
        diff_T_mean100(noise_level,iterative_index)=norm(T_rot_average-T_real);
        %         diff_T_kneip(noise_level,iterative_index)=norm(T_Kneip-T_real);
        %         diff_T_Fischler(noise_level,iterative_index)=norm(T_fishler-T_real);
        %         diff_T_QuanLan(noise_level,iterative_index)=norm(T_LongQun-T_real);
        %         diff_T_QuanLan_mp(noise_level,iterative_index)=norm(T_LongQun_mp-T_real);
        
        
        %         diff_T_3p(noise_level,iterative_index)=norm(tp3-Rp3'*T_real');
        %         diff_T_4p(noise_level,iterative_index)=norm(tp4-Rp4'*T_real');
        %         diff_T_epnp(noise_level,iterative_index)=norm(Tp-T_real');
        %         diff_T_Paul(noise_level,iterative_index)=norm(TPaul-RPaul'*T_real');
        disp(num2str(iterative_index));
    end
end
figure(1)
plot(diff_Rot_Paul(10,:),'.');
axis([0,100,0,2])
title('Fiore');
ylabel('Rotation error')
% figure(2)
% plot(diff_Rot_kneip(10,:),'.');
% axis([0,100,0,2])
% title('Kneip');
% ylabel('Rotation error')
% % hold on%         diff_Rot_3p(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Rp3-eye(3))));
% %         diff_Rot_4p(noise_level,iterative_index)=sum(sum(abs(Rot_origin*Rp4-eye(3))));
% figure(3)
% plot(diff_Rot_Fischler(10,:),'.');
% axis([0,100,0,2])
% title('Fischler & Bolles');
% ylabel('Rotation error')
% figure(4)
% % hold on
% plot(diff_Rot_QuanLan(10,:),'.');
% axis([0,100,0,2])
% title('Quan & Lan');
% ylabel('Rotation error')
% % hold on
% figure(5)
% plot(diff_Rot_QuanLan_mp(10,:),'.');
% axis([0,100,0,2])
% title('Quan & Lan multi-precision');
% ylabel('Rotation error')
figure(6)
plot(diff_Rot_epnp(10,:),'.');
axis([0,100,0,2])
title('Lepetit');
ylabel('Rotation error')
figure(7)
plot(diff_Rot_mean50(10,:),'.');
axis([0,100,0,2])
title('Rotation Averaging 50');
ylabel('Rotation error')
figure(8)
plot(diff_Rot_mean100(10,:),'.');
axis([0,100,0,2])
title('Rotation Averaging 100');
ylabel('Rotation error')
% hold on

figure(11)
plot(diff_T_Paul(10,:),'.');
axis([0,100,0,20])
title('Fiore');
ylabel('Translation error')
% figure(12)
% plot(diff_T_kneip(10,:),'.');
% axis([0,100,0,10])
% title('Kneip');
% ylabel('Translation error')
% % hold on
% figure(13)
% plot(diff_T_Fischler(10,:),'.');
% axis([0,100,0,10])
% title('Fischler & Bolles');
% ylabel('Translation error')
% figure(14)
% % hold on
% plot(diff_T_QuanLan(10,:),'.');
% axis([0,100,0,10])
% title('Quan & Lan');
% ylabel('Translation error')
% % hold on
% figure(15)
% plot(diff_T_QuanLan_mp(10,:),'.');
% axis([0,100,0,10])
% title('Quan & Lan multi-precision');
% ylabel('Translation error')
% pause(0.1)
figure(16)
plot(diff_T_epnp(10,:),'.');
axis([0,100,0,20])
title('Lepetit');
ylabel('Translation error')
figure(17)
plot(diff_T_mean50(10,:),'.');
axis([0,100,0,20])
title('Rotation Averaging 50');
ylabel('Translation error')
figure(18)
plot(diff_T_mean100(10,:),'.');
axis([0,100,0,20])
title('Rotation Averaging 100');
ylabel('Translation error')
pause(0.1)

figure(30)
plot(mean(diff_T_Paul,2),'*')
hold on
plot(mean(diff_T_epnp,2),'.')
% hold on
% plot(mean(diff_T_epnp,2),'.')
hold on
plot(mean(diff_T_mean50,2),'o')
hold on
plot(mean(diff_T_mean100,2),'x')
% hold on
% plot(mean(diff_T_Kneip,2),'o')
ylabel('Mean Translation error')
legend('Fiore','Lepetit','RA50','RA100')
axis([1,10,0,5])

figure(31)
plot(median(diff_T_Paul,2),'*')
hold on
plot(median(diff_T_epnp,2),'.')
% hold on
% plot(mean(diff_T_epnp,2),'.')
hold on
plot(median(diff_T_mean50,2),'o')
hold on
plot(median(diff_T_mean100,2),'x')
% hold on
% plot(mean(diff_T_Kneip,2),'o')
ylabel('Median Translation error')
legend('Fiore','Lepetit','RA50','RA100')
axis([1,10,0,5])

figure(32)
plot(mean(diff_Rot_Paul,2),'*')
hold on
plot(mean(diff_Rot_epnp,2),'.')
% hold on
% plot(mean(diff_R_epnp,2),'.')
hold on
plot(mean(diff_Rot_mean50,2),'o')
hold on
plot(mean(diff_Rot_mean100,2),'x')
% hold on
% plot(mean(diff_T_Kneip,2),'o')
ylabel('Mean Rotation error')
legend('Fiore','Lepetit','RA50','RA100')
axis([1,10,0,0.5])

figure(33)
plot(median(diff_Rot_Paul,2),'*')
hold on
plot(median(diff_Rot_epnp,2),'.')
% hold on
% plot(mean(diff_T_epnp,2),'.')
hold on
plot(median(diff_Rot_mean50,2),'o')
hold on
plot(median(diff_Rot_mean100,2),'x')
% hold on
% plot(mean(diff_T_Kneip,2),'o')
ylabel('Median Rotation error')
legend('Fiore','Lepetit','RA50','RA100')
axis([1,10,0,0.5])

save data_PNP_C2_100 diff_T_mean50 diff_T_mean100 diff_Rot_mean50 diff_Rot_mean100 diff_T_epnp diff_T_Paul diff_Rot_epnp diff_Rot_Paul













