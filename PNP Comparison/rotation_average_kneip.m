function [R_rot_average,T_rot_average]=rotation_average_kneip(poses1,poses2,poses3,poses4)
R_rot_average=eye(3);
T_rot_average=[0,0,0]';

R11=poses1(:,2:4);
R12=poses1(:,6:8);
R13=poses1(:,10:12);
R14=poses1(:,14:16);

R21=poses2(:,2:4);
R22=poses2(:,6:8);
R23=poses2(:,10:12);
R24=poses2(:,14:16);

R31=poses3(:,2:4);
R32=poses3(:,6:8);
R33=poses3(:,10:12);
R34=poses3(:,14:16);

R41=poses4(:,2:4);
R42=poses4(:,6:8);
R43=poses4(:,10:12);
R44=poses4(:,14:16);

T(1,1,:)=-R11'*poses1(:,1);
T(1,2,:)=-R12'*poses1(:,5);
T(1,3,:)=-R13'*poses1(:,9);
T(1,4,:)=-R14'*poses1(:,13);

T(2,1,:)=-R21'*poses2(:,1);
T(2,2,:)=-R22'*poses2(:,5);
T(2,3,:)=-R23'*poses2(:,9);
T(2,4,:)=-R24'*poses2(:,13);

T(3,1,:)=-R31'*poses3(:,1);
T(3,2,:)=-R32'*poses3(:,5);
T(3,3,:)=-R33'*poses3(:,9);
T(3,4,:)=-R34'*poses3(:,13);

T(4,1,:)=-R41'*poses4(:,1);
T(4,2,:)=-R42'*poses4(:,5);
T(4,3,:)=-R43'*poses4(:,9);
T(4,4,:)=-R44'*poses4(:,13);

T11=-R11'*poses1(:,1);
T12=-R12'*poses1(:,5);
T13=-R13'*poses1(:,9);
T14=-R14'*poses1(:,13);

T21=-R21'*poses2(:,1);
T22=-R22'*poses2(:,5);
T23=-R23'*poses2(:,9);
T24=-R24'*poses2(:,13);

T31=-R31'*poses3(:,1);
T32=-R32'*poses3(:,5);
T33=-R33'*poses3(:,9);
T34=-R34'*poses3(:,13);

T41=-R41'*poses4(:,1);
T42=-R42'*poses4(:,5);
T43=-R43'*poses4(:,9);
T44=-R44'*poses4(:,13);

[theta11, vector_y11]=rodrigues_rot2vetor(R11) ;
[theta12, vector_y12]=rodrigues_rot2vetor(R12); 
[theta13, vector_y13]=rodrigues_rot2vetor(R13) ;
[theta14, vector_y14]=rodrigues_rot2vetor(R14); 

[theta21, vector_y21]=rodrigues_rot2vetor(R21) ;
[theta22, vector_y22]=rodrigues_rot2vetor(R22); 
[theta23, vector_y23]=rodrigues_rot2vetor(R23) ;
[theta24, vector_y24]=rodrigues_rot2vetor(R24);

[theta31, vector_y31]=rodrigues_rot2vetor(R31) ;
[theta32, vector_y32]=rodrigues_rot2vetor(R32); 
[theta33, vector_y33]=rodrigues_rot2vetor(R33) ;
[theta34, vector_y34]=rodrigues_rot2vetor(R34); 

[theta41, vector_y41]=rodrigues_rot2vetor(R41) ;
[theta42, vector_y42]=rodrigues_rot2vetor(R42); 
[theta43, vector_y43]=rodrigues_rot2vetor(R43) ;
[theta44, vector_y44]=rodrigues_rot2vetor(R44); 
p(1,1,:)=theta11* vector_y11;
p(1,2,:)=theta12* vector_y12;
p(1,3,:)=theta13* vector_y13;
p(1,4,:)=theta14* vector_y14;

p(2,1,:)=theta21* vector_y21;
p(2,2,:)=theta22* vector_y22;
p(2,3,:)=theta23* vector_y23;
p(2,4,:)=theta24* vector_y24;

p(3,1,:)=theta31* vector_y31;
p(3,2,:)=theta32* vector_y32;
p(3,3,:)=theta33* vector_y33;
p(3,4,:)=theta34* vector_y34;


p(4,1,:)=theta41* vector_y41;
p(4,2,:)=theta42* vector_y42;
p(4,3,:)=theta43* vector_y43;
p(4,4,:)=theta44* vector_y44;

p11=theta11* vector_y11;
p12=theta12* vector_y12;
p13=theta13* vector_y13;
p14=theta14* vector_y14;

p21=theta21* vector_y21;
p22=theta22* vector_y22;
p23=theta23* vector_y23;
p24=theta24* vector_y24;

p31=theta31* vector_y31;
p32=theta32* vector_y32;
p33=theta33* vector_y33;
p34=theta34* vector_y34;


p41=theta41* vector_y41;
p42=theta42* vector_y42;
p43=theta43* vector_y43;
p44=theta44* vector_y44;

min_distance=100000;
for solution_index_group1=1:4
    for solution_index_group2=1:4
        for solution_index_group3=1:4
            for solution_index_group4=1:4
                Ta=strcat('T1',num2str(solution_index_group1));
                distance_set_T=...
                    norm([T(1,solution_index_group1,1)-T(2,solution_index_group2,1),T(1,solution_index_group1,2)-T(2,solution_index_group2,2),T(1,solution_index_group1,3)-T(2,solution_index_group2,3)])+...
                    norm([T(1,solution_index_group1,1)-T(3,solution_index_group3,1),T(1,solution_index_group1,2)-T(3,solution_index_group3,2),T(1,solution_index_group1,3)-T(3,solution_index_group3,3)])+...
                    norm([T(1,solution_index_group1,1)-T(4,solution_index_group4,1),T(1,solution_index_group1,2)-T(4,solution_index_group4,2),T(1,solution_index_group1,3)-T(4,solution_index_group4,3)])+...
                    norm([T(2,solution_index_group2,1)-T(3,solution_index_group3,1),T(2,solution_index_group2,2)-T(3,solution_index_group3,2),T(2,solution_index_group2,3)-T(3,solution_index_group3,3)])+...
                    norm([T(2,solution_index_group2,1)-T(4,solution_index_group4,1),T(2,solution_index_group2,2)-T(4,solution_index_group4,2),T(2,solution_index_group2,3)-T(4,solution_index_group4,3)])+...
                    norm([T(3,solution_index_group3,1)-T(4,solution_index_group4,1),T(3,solution_index_group3,2)-T(4,solution_index_group4,2),T(3,solution_index_group3,3)-T(4,solution_index_group4,3)]);
               distance_set_p=...
                    norm([p(1,solution_index_group1,1)-p(2,solution_index_group2,1),p(1,solution_index_group1,2)-p(2,solution_index_group2,2),p(1,solution_index_group1,3)-p(2,solution_index_group2,3)])+...
                    norm([p(1,solution_index_group1,1)-p(3,solution_index_group3,1),p(1,solution_index_group1,2)-p(3,solution_index_group3,2),p(1,solution_index_group1,3)-p(3,solution_index_group3,3)])+...
                    norm([p(1,solution_index_group1,1)-p(4,solution_index_group4,1),p(1,solution_index_group1,2)-p(4,solution_index_group4,2),p(1,solution_index_group1,3)-p(4,solution_index_group4,3)])+...
                    norm([p(2,solution_index_group2,1)-p(3,solution_index_group3,1),p(2,solution_index_group2,2)-p(3,solution_index_group3,2),p(2,solution_index_group2,3)-p(3,solution_index_group3,3)])+...
                    norm([p(2,solution_index_group2,1)-p(4,solution_index_group4,1),p(2,solution_index_group2,2)-p(4,solution_index_group4,2),p(2,solution_index_group2,3)-p(4,solution_index_group4,3)])+...
                    norm([p(3,solution_index_group3,1)-p(4,solution_index_group4,1),p(3,solution_index_group3,2)-p(4,solution_index_group4,2),p(3,solution_index_group3,3)-p(4,solution_index_group4,3)]);
%                 isall_valide=length(find([solution_set_1,solution_set_2,solution_set_3,solution_set_4]==0));
                distance_this_set=mean(distance_set_T)*mean(distance_set_p);
                if distance_this_set<min_distance
                    min_distance=distance_this_set;
                    best_set_index=[solution_index_group1,solution_index_group2,solution_index_group3,solution_index_group4];
                end
            end
        end
    end
end
T_rot_average1=(T(1,best_set_index(1),:)+T(2,best_set_index(2),:)+T(3,best_set_index(3),:)+T(4,best_set_index(4),:))/4;
T_rot_average=[T_rot_average1(1,1,1),T_rot_average1(1,1,2),T_rot_average1(1,1,3)]';
distance12v=p(1,best_set_index(1),:)-p(2,best_set_index(2),:);
distance12=norm([distance12v(1,1,1),distance12v(1,1,2),distance12v(1,1,3)]);
distance13v=p(1,best_set_index(1),:)-p(3,best_set_index(3),:);
distance13=norm([distance13v(1,1,1),distance13v(1,1,2),distance13v(1,1,3)]);
distance14v=p(1,best_set_index(1),:)-p(4,best_set_index(4),:);
distance14=norm([distance14v(1,1,1),distance14v(1,1,2),distance14v(1,1,3)]);
distance23v=p(2,best_set_index(2),:)-p(3,best_set_index(3),:);
distance23=norm([distance23v(1,1,1),distance23v(1,1,2),distance23v(1,1,3)]);
distance24v=p(2,best_set_index(2),:)-p(4,best_set_index(4),:);
distance24=norm([distance24v(1,1,1),distance24v(1,1,2),distance24v(1,1,3)]);
distance34v=p(3,best_set_index(3),:)-p(4,best_set_index(4),:);
distance34=norm([distance34v(1,1,1),distance34v(1,1,2),distance34v(1,1,3)]);

distance1=distance12+distance13+distance14;
distance2=distance12+distance23+distance24;
distance3=distance13+distance23+distance34;
distance4=distance14+distance24+distance34;
[max_distance,max_index]=max([distance1,distance2,distance3,distance4]);
reserve_index=setdiff([1,2,3,4],[max_index]);
p_average1=(p(reserve_index(1),best_set_index(reserve_index(1)),:)+p(reserve_index(2),best_set_index(reserve_index(2)),:)+p(reserve_index(3),best_set_index(reserve_index(3)),:))/3;
% p_average1=(p(1,best_set_index(1),:)+p(2,best_set_index(2),:)+p(3,best_set_index(3),:)+p(4,best_set_index(4),:))/4;
p_average=[p_average1(1,1,1),p_average1(1,1,2),p_average1(1,1,3)]';
theta_average=norm(p_average);
vector_average=p_average/norm(p_average);
R_rot_average=rodrigues_vetor2rot(theta_average,vector_average)' ;
%  R_Kneip=Rp3'
%  T_Kneip=-Rp3'*tp3;
% distance


