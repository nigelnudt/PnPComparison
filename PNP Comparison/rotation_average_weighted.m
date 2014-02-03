function [R_rot_average,T_rot_average]=rotation_average_weighted(R_rot_average,T_rot_average,pose,weight)
[theta_avg, vector_avg]=rodrigues_rot2vetor(R_rot_average) ;
p_avg=theta_avg* vector_avg;
R1=pose(:,2:4)';
R2=pose(:,6:8)';
R3=pose(:,10:12)';
R4=pose(:,14:16)';
T1=-R1*pose(:,1);
T2=-R2*pose(:,5);
T3=-R3*pose(:,9);
T4=-R4*pose(:,13);
[theta1, vector_y1]=rodrigues_rot2vetor(R1) ;
[theta2, vector_y2]=rodrigues_rot2vetor(R2);
[theta3, vector_y3]=rodrigues_rot2vetor(R3) ;
[theta4, vector_y4]=rodrigues_rot2vetor(R4);
p1=theta1* vector_y1;
p2=theta2* vector_y2;
p3=theta3* vector_y3;
p4=theta4* vector_y4;
for index=1:4
    r_name=strcat('p',num2str(index));
    distance_r(index)=norm(p_avg-eval(r_name));
end
for index=5:8
    r_name=strcat('p',num2str(index-4));
    distance_r(index)=norm(-p_avg-eval(r_name));
end
[min_r_distance,index_select]=min(distance_r);
if index_select>4
    index_select=index_select-4;
    turn_over=-1;
else
    index_select=index_select;
    turn_over=1;
end
p_this=eval(strcat('p',num2str(index_select)));
if min_r_distance<5
    p_avg=turn_over*(1-weight)*p_avg+weight*p_this;
    theta_average=norm(p_avg);
    vector_average=p_avg/norm(p_avg);
    R_rot_average=rodrigues_vetor2rot(theta_average,vector_average)' ;
    T_this=eval(strcat('T',num2str(index_select)));
    T_rot_average=(1-weight)*T_rot_average+weight*T_this;
else
    R_rot_average=R_rot_average;
    T_rot_average=T_rot_average;
end