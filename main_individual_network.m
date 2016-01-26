g = genpath('../'); addpath(g);
curfold=pwd;
close all;
% Initialize random number generator
rng(37695,'twister'); 
% Initialize parameters
n=50; % network size
pH=0; % prob. of losing resistance
pV=0; % prob. of losing infectivity
pr=1; % prob. of developing resistance (assumed to be 1)
p0i=1; % prob. of developing infectivity (assumed to be 1)
istart=1; inet=1; % size of and type of initial network
%inet= 1. for one-to-one or diagonal network initialization
%2. for perfectly modular initial network with two modules
%3. for perfectly nested initial network 
%4. for network that minimizes the sum of modularity and nestedness at size 4x4 and connectance 0.5
%5. Modified Erdos-Renyi random network with all diagonal elements equal to 1
%and average overall connectance specified by pnet
pnet=0.5; % connection probability for modified ER random network
% Initialize network
if inet==1 % one-to-one IC
    A0=eye(istart);
elseif inet==2 % Modular IC
    A0=zeros(istart);
    imid=ceil(istart/2);
    A0(1:imid,1:imid)=1; A0(imid+1:istart,imid+1:istart)=1;
elseif inet==3 % Nested IC
    A0=triu(ones(istart));
elseif inet==4 % Minimized IC
    istart=4;
    A0=zeros(4);A0(1,3:4)=1;A0(2,1:2)=1;A0(3,2)=1;A0(3,4)=1;A0(4,1)=1;A0(4,3)=1;
elseif inet==5 % Modified ER random IC
    peff=(pnet*istart-1)/(istart-1);
    A0=rand(istart); A0=A0<peff;
    A0(1:istart+1:end)=ones(1,istart); 
end
% Define variables to be computed
global A kh kv % Adjacency matrix of network, host degree and virus degree
tic
i=istart;
iv=istart;
A=zeros(n);
A(1:i,1:i)=A0;
A(1,1)=1; 
kh=sum(A,2);
kv=sum(A,1);
% Network growth
while i<n
    % Determine which nodes to replicate
    khred=kh(1:i); 
    kvred=kv(1:iv);
    [ihmin ihmax]=findminmax(khred,i);
    [ivmin ivmax]=findminmax(kvred,iv);
    % Add new host and rewire
        i=i+1;
        ih=i;
    for j=1:iv
        Aold=A(ihmin,j);
        rewire(ih,j,Aold,1.0-pr,pH);
    end
    % Add new virus and rewire
    iv=iv+1;
    for j=1:i-1
        Aold=A(j,ivmax);
        rewire(j,iv,Aold,1.0-pV,p0i);
    end
    Aold=A(ih,ivmax);
    rewire(i,iv,Aold,1.0-pV,p0i);
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