%% Output of parameters and nestedness/modularity
% For conflicting attachment model
inet=1;
dlmwrite([curfold '/n.txt'],para,' ')
dlmwrite([curfold '/nestedness_vs_n_conflict.txt'],[para nest(:,inet) nest_std(:,inet)],' ')
dlmwrite([curfold '/modularity_Qb_vs_n_conflict.txt'],[para modu(:,inet) modu_std(:,inet)],' ')
dlmwrite([curfold '/modularity_Qr_vs_n_conflict.txt'],[para ierat(:,inet) ierat_std(:,inet)],' ')
dlmwrite([curfold '/modularity_comm_no_vs_n_conflict.txt'],[para ncom(:,inet) ncom_std(:,inet)],' ')
dlmwrite([curfold '/module_nestedness_zscore_vs_n_conflict.txt'],[para znest znest_std],' ')
dlmwrite([curfold '/module_zscore_usable_ratio_vs_n_conflict.txt'],[para nmod_rat],' ')
% For random network
inet=2;
dlmwrite([curfold '/n.txt'],para,' ')
dlmwrite([curfold '/nestedness_vs_n_random.txt'],[para nest(:,inet) nest_std(:,inet)],' ')
dlmwrite([curfold '/modularity_Qb_vs_n_random.txt'],[para modu(:,inet) modu_std(:,inet)],' ')
dlmwrite([curfold '/modularity_Qr_vs_n_random.txt'],[para ierat(:,inet) ierat_std(:,inet)],' ')
dlmwrite([curfold '/modularity_comm_no_vs_n_random.txt'],[para ncom(:,inet) ncom_std(:,inet)],' ')