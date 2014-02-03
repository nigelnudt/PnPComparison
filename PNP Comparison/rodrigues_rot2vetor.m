function [theta, vector_y]=rodrigues_rot2vetor(Rot_origin) 
y_asym_matrix=(Rot_origin-Rot_origin')/2;
vector_y=[y_asym_matrix(3,2),-y_asym_matrix(3,1),y_asym_matrix(2,1)];
theta=asin(norm(vector_y));
vector_y=vector_y/norm(vector_y);
