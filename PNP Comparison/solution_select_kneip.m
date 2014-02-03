function [R_Kneip,T_Kneip]=solution_select_kneip(poses,X3d_world,x2d_h_normlized)
R_Kneip=eye(3);
T_Kneip=zeros(3,1);
verify_error=100;
for index=1:4
    Rp3=poses(:,4*index-2:4*index);
    tp3=poses(:,4*index-3);
    p4_camera=Rp3'*(X3d_world'-tp3);
    p4_error=norm([p4_camera(1)/p4_camera(3),p4_camera(2)/p4_camera(3)]-x2d_h_normlized(1:2));
    if p4_error<verify_error
        R_Kneip=Rp3';
        T_Kneip=-Rp3'*tp3;
        verify_error=p4_error;
    end
end