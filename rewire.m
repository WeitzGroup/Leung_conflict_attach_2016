function []=rewire(i,j,Aold,p11,p01)
% Determine the connectivity of a newly added node
global A kh kv
randum=rand;
if Aold==1
    det=randum<p11;
else
    det=randum<p01;
end
A(i,j)=det; 
kh(i)=kh(i)+det;
kv(j)=kv(j)+det;