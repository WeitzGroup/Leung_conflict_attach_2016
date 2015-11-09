function [ikmin ikmax]=findminmax(k,i)
% Find nodes with the minimum and maximum degree
[sk,isk]=sort(k);
kmax=sk(i);
kmin=sk(1);
iskmax=k==kmax;
iskmin=k==kmin;
ikmax=isk(i+1-unidrnd(sum(iskmax)));
ikmin=isk(unidrnd(sum(iskmin)));