g = genpath('../'); addpath(g);
curfold=pwd;
close all;
% Initialize parameters
net_type=1; % 1 for conflicting attachment model, 0 for random network
pran=0.5; % Connection prob. of random network (if chosen)
nmin=10; nmax=100; nstep=1; % network size sampled from nmin to nmax with interval nstep
n=nmin; nt=int64(1+(nmax-nmin)/nstep);
pH=0; % prob. of losing resistance
pV=0; % prob. of losing infectivity
pr=1; % prob. of developing resistance (assumed to be 1)
p0i=1; % prob. of developing infectivity (assumed to be 1)
kth=0; % Virus degree below which to remove at fixed host no. stage
nfix=0; % no. of steps to run with fixed host no. (0 means no host/virus extinction stage)
para=transpose(nmin:nstep:nmax); % matrix containing the parameter n tested
ient=10000.*(10./para).^2; ient(ient<1)=1; % no. of networks in ensemble average, rounded to nearest +ve integer
ient=int64(ient); dent=double(ient);
% Define variables to be computed
global A kh kv % Adjacency matrix of network, host degree and virus degree
nest=zeros(nt,1); % nestedness
modu=zeros(nt,1); % modularity
ncom=zeros(nt,1); % no. of modules
ierat=zeros(nt,1); % ratio of interior to exterior edges
asp=zeros(nt,1); % aspect ratio
% Initialize random number generator
%rng(15614,'twister'); 
%rng(17589,'twister'); 
rng(37695,'twister'); 
tic
% Loop through different network sizes
for in=1:nt
        for ien=1:ient(in)
% Conflicting attachment model
if net_type==1
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
% Random network model
elseif net_type==0
i=n; iv=n;
A=double(rand(n)<pran);
end
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
modu(in)=modu(in)+bp.community.Qb;
ierat(in)=ierat(in)+bp.community.Qr;
ncom(in)=ncom(in)+bp.community.N;
asp(in)=asp(in)+i/iv;
% Nestedness
bp.nestedness.Detect();
nest(in)=nest(in)+bp.nestedness.N;
%% Network plots
%run('netplot')
        end
        % Compute ensemble averages
        dentemp=dent(in);
        modu(in)=modu(in)/dentemp;
        ierat(in)=ierat(in)/dentemp;
        ncom(in)=ncom(in)/dentemp;
        nest(in)=nest(in)/dentemp;
        asp(in)=asp(in)/dentemp;
% Change network size n
n=n+nstep;
end
toc
% Output of the modularity/nestedness information
run('sizeout')
% Plots of modularity/nestedness vs network size
run('sizeplot')