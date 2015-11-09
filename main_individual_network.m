g = genpath('../'); addpath(g);
curfold=pwd;
close all;
% Initialize parameters
n=50; % network size
pH=0; % prob. of losing resistance
pV=0; % prob. of losing infectivity
pr=1; % prob. of developing resistance (assumed to be 1)
p0i=1; % prob. of developing infectivity (assumed to be 1)
kth=0; % Virus degree below which to remove at fixed host no. stage
nfix=0; % no. of steps to run with fixed host no. (0 means no host/virus extinction stage)
% Define variables to be computed
global A kh kv % Adjacency matrix of network, host degree and virus degree
% Initialize random number generator
%rng(15614,'twister'); 
%rng(17589,'twister'); 
rng(37695,'twister'); 
tic
A=zeros(n+nfix);
i=1;
ik=0;
iv=1;
ifix=0;
A(1,1)=1; 
kh=sum(A,2);
kv=sum(A,1);
% Network growth
while ifix<nfix || i<n || ik==1
    % Determine which nodes to replicate
    khred=kh(1:i); 
    kvred=kv(1:iv);
    [ihmin ihmax]=findminmax(khred,i);
    [ivmin ivmax]=findminmax(kvred,iv);
    % Add new host and rewire
    if i<n || ik==1
        i=i+1;
        ih=i;
    for j=1:iv
        Aold=A(ihmin,j);
        rewire(ih,j,Aold,1.0-pr,pH);
    end
    % Add new virus and rewire
    if ik==0 %|| kh(ih)==0
    iv=iv+1;
    for j=1:i-1
        Aold=A(j,ivmax);
        rewire(j,iv,Aold,1.0-pV,p0i);
    end
    Aold=A(ih,ivmax);
    rewire(i,iv,Aold,1.0-pV,p0i);
    end
    ik=0;
    elseif nfix~=0
        % Extinction of host and virus
        ifix=ifix+1;
        vrange=1:iv;
        ivb=A(ihmax,vrange)==1;
        kv(vrange)=kv(vrange)-ivb;
        A(ihmax,:)=[];
        kh(ihmax)=[];
        i=i-1;
        ivd=kv(vrange)<=kth;
        A(:,ivd)=[];
        kv(ivd)=[];
        iv=iv-sum(ivd);
        ik=1;
    end
end
A=A(1:i,1:iv);
%% Output connectivity matrix
%dlmwrite([curfold '/network.txt'],A,' ')
%% Bipartite matrix analysis
bp=Bipartite(A);
%bp.row_labels=row_labels;
%bp.col_labels=col_labels;
%bp.printer.PrintGeneralProperties();
% Community structure
bp.community = AdaptiveBrim(bp.matrix);
bp.community.Detect();
modu=bp.community.Qb;
ierat=bp.community.Qr;
ncom=bp.community.N;
asp=i/iv;
% Nestedness
bp.nestedness.Detect();
nest=bp.nestedness.N;
%% Network plots
run('netplot')
toc