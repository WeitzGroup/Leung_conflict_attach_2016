%% Output of parameters and nestedness/modularity
dlmwrite([curfold '/pl.txt'],para(:,:,1),' ')
dlmwrite([curfold '/pf.txt'],para(:,:,2),' ')
dlmwrite([curfold '/nestedness.txt'],nest,' ')
dlmwrite([curfold '/modularity_Qb.txt'],modu,' ')
dlmwrite([curfold '/modularity_Qr.txt'],ierat,' ')
dlmwrite([curfold '/modularity_comm_no.txt'],ncom,' ')