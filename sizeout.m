%% Output of parameters and nestedness/modularity
dlmwrite([curfold '/n.txt'],para,' ')
dlmwrite([curfold '/nestedness_vs_n.txt'],[para nest],' ')
dlmwrite([curfold '/modularity_Qb_vs_n.txt'],[para modu],' ')
dlmwrite([curfold '/modularity_Qr_vs_n.txt'],[para ierat],' ')
dlmwrite([curfold '/modularity_comm_no_vs_n.txt'],[para ncom],' ')
dlmwrite([curfold '/aspect_ratio_vs_n.txt'],[para asp],' ')