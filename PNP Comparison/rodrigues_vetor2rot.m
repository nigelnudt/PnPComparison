function [Rot_recover]=rodrigues_vetor2rot(theta, vector_y)      
vector_v=vector_y/norm(vector_y);
v_asym_matrix=[0,-vector_v(3),vector_v(2);vector_v(3),0,-vector_v(1);-vector_v(2),vector_v(1),0];
Rot_recover=eye(3)+sin(theta)*v_asym_matrix+(1-cos(theta))*v_asym_matrix*v_asym_matrix;