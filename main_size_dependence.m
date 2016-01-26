g = genpath('../'); addpath(g);
curfold=pwd;
close all;
% Initialize random number generator
rng(37695,'twister'); 
% Initialize parameters
nmin=10; nmax=100; nstep=10; % network size sampled from nmin to nmax with interval nstep
n=nmin; nt=int64(1+(nmax-nmin)/nstep);
pH=0.8; % prob. of losing resistance
pV=0.2; % prob. of losing infectivity
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
end
para=transpose(nmin:nstep:nmax); % matrix containing the parameter n tested
ient=10000.*(10./para).^2; ient(ient<1)=1; % no. of networks in ensemble average, rounded to nearest +ve integer
ient=int64(ient); net_types=2;
% Define variables to be computed
global A kh kv % Adjacency matrix of network, host degree and virus degree
nest=zeros(nt,net_types); nest_std=zeros(nt,net_types); % nestedness and its standard deviation
modu=zeros(nt,net_types); modu_std=zeros(nt,net_types); % modularity and its std
ncom=zeros(nt,net_types); ncom_std=zeros(nt,net_types); % no. of modules and its std
ierat=zeros(nt,net_types); ierat_std=zeros(nt,net_types); % ratio of interior to exterior edges and its std
znest=zeros(nt,1); znest_std=zeros(nt,1); % z-score of nestedness within modules and its std
nmod_rat=zeros(nt,1); % Ratio of usable modules (not fully connected) when calculating the z-scores
tic
% Loop through different network sizes
for in=1:nt
nest_en=zeros(ient(in),net_types);
modu_en=zeros(ient(in),net_types);
ncom_en=zeros(ient(in),net_types);
ierat_en=zeros(ient(in),net_types);
nmod_allo=5*ient(in);
znest_en=zeros(nmod_allo,1);
nmod=0; totmod=0;
        for ien=1:ient(in)
% Conflicting attachment model
if inet==5 % Modified ER random IC
    peff=(pnet*istart-1)/(istart-1);
    A0=rand(istart); A0=A0<peff;
    A0(1:istart+1:end)=ones(1,istart); 
end
A=zeros(n);
i=istart;
iv=istart;
A(1:i,1:i)=A0; 
kh=sum(A,2);
kv=sum(A,1);
% Network growth
while i<n
    % Determine which nodes to replicate
    khred=kh(1:i); 
    kvred=kv(1:iv);
    [ihmin, ihmax]=findminmax(khred,i);
    [ivmin, ivmax]=findminmax(kvred,iv);
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
% Random network model with the same connectance
pran=sum(sum(A))/(i*iv);
R=double(rand(i,iv)<pran);
AR=zeros(i,iv,net_types);
AR(:,:,1)=A; AR(:,:,2)=R;
%% Bipartite matrix analysis
for inet=1:net_types
bp=Bipartite(AR(:,:,inet));
% Community structure
bp.community = AdaptiveBrim(bp.matrix);
bp.community.Detect();
modu_en(ien,inet)=bp.community.Qb;
ierat_en(ien,inet)=bp.community.Qr;
ncom_en(ien,inet)=bp.community.N;
% Nestedness
bp.nestedness.Detect();
nest_en(ien,inet)=bp.nestedness.N;
% Statistics within modules (multi-scale structure)
if inet==1
bp.internal_statistics.TestInternalModules(100,@NullModels.EQUIPROBABLE);
zn=bp.internal_statistics.meta_statistics.N_values.zscore;
totmod=totmod+bp.community.N;
for imod=1:bp.community.N
if ~isnan(zn(imod))
nmod=nmod+1;
znest_en(nmod)=zn(imod);
end
end
end
end
        end
        % Compute ensemble averages
        nmod_range=1:nmod;
        znest_en=znest_en(nmod_range);
        nmod_rat(in,:)=nmod/totmod;
        modu(in,:)=mean(modu_en);
        ierat(in,:)=mean(ierat_en);
        ncom(in,:)=mean(ncom_en);
        nest(in,:)=mean(nest_en);
        znest(in)=mean(znest_en);
        modu_std(in,:)=std(modu_en);
        ierat_std(in,:)=std(ierat_en);
        ncom_std(in,:)=std(ncom_en);
        nest_std(in,:)=std(nest_en);
        znest_std(in)=std(znest_en);
% Change network size n
n=n+nstep;
end
toc
% Output of the modularity/nestedness information
run('sizeout')
% Plots of modularity/nestedness vs network size
run('sizeplot')