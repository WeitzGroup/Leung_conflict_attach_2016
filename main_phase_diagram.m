g = genpath('../'); addpath(g);
curfold=pwd;
close all;
% Initialize parameters
n=50; % network size
pH=0; % prob. of losing resistance (start sampling from 0)
pV=0; % prob. of losing infectivity (start sampling from 0)
pr=1; % prob. of developing resistance (assumed to be 1)
p0i=1; % prob. of developing infectivity (assumed to be 1)
kth=0; % Virus degree below which to remove at fixed host no. stage
nfix=0; % no. of steps to run with fixed host no. (0 means no host/virus extinction stage)
ient=20; % no. of networks in ensemble average
pstep=0.01; % Increments of parameters used in phase diagram
nt=int64(1+1/pstep); % phase diagram has nt*nt no. of points
para=repmat(0,[nt nt 2]); % matrix containing the parameters pH and pV
% Define variables to be computed
global A kh kv % Adjacency matrix of network, host degree and virus degree
nest=zeros(nt); % nestedness
modu=zeros(nt); % modularity
ncom=zeros(nt); % no. of modules
ierat=zeros(nt); % ratio of interior to exterior edges
asp=zeros(nt); % aspect ratio
% Initialize random number generator
%rng(15614,'twister'); 
%rng(17589,'twister'); 
rng(37695,'twister'); 
tic
% Loop through pH
for ipH=1:nt
    pV=0;
    % Loop through pV
    for ipV=1:nt
        para(ipH,ipV,1)=pH;
        para(ipH,ipV,2)=pV;
        for ien=1:ient
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
modu(ipH,ipV)=modu(ipH,ipV)+bp.community.Qb;
ierat(ipH,ipV)=ierat(ipH,ipV)+bp.community.Qr;
ncom(ipH,ipV)=ncom(ipH,ipV)+bp.community.N;
asp(ipH,ipV)=asp(ipH,ipV)+i/iv;
% Nestedness
bp.nestedness.Detect();
nest(ipH,ipV)=nest(ipH,ipV)+bp.nestedness.N;
%% Network plots
%run('netplot')
        end
        % Compute ensemble averages
        modu(ipH,ipV)=modu(ipH,ipV)/ient;
        ierat(ipH,ipV)=ierat(ipH,ipV)/ient;
        ncom(ipH,ipV)=ncom(ipH,ipV)/ient;
        nest(ipH,ipV)=nest(ipH,ipV)/ient;
        asp(ipH,ipV)=asp(ipH,ipV)/ient;
    % Change pV
    pV=pV+pstep;
    end
% Change pH
pH=pH+pstep;
end
toc
% Output of the modularity/nestedness information
run('phaseout')
% Color map plots of the modularity/nestedness
run('phaseplot')