%% Plotting the dependence of modularity/nestedness on network size
% Conflicting attachment model
inet=1;
% Modularity
figure(1);
errorbar(para,modu(:,inet),modu_std(:,inet))
xlabel('n'); ylabel('Modularity Qb');
savefig('modularity_Qb_vs_n_conflict')
% Nestedness
figure(2); 
errorbar(para,nest(:,inet),nest_std(:,inet))
xlabel('n'); ylabel('Nestedness');
savefig('nestedness_vs_n_conflict')
% Interia/exteria edge ratio
figure(3);
errorbar(para,ierat(:,inet),ierat_std(:,inet))
xlabel('n'); ylabel('Interia to exteria links ratio');
savefig('modularity_Qr_vs_n_conflict')
% Module nestedness z-score
figure(4); 
errorbar(para,znest,znest_std)
xlabel('n'); ylabel('Module nestedness z-score');
savefig('module_nestedness_zscore_vs_n_conflict')
%Random network
inet=2;
% Modularity
figure(5);
errorbar(para,modu(:,inet),modu_std(:,inet))
xlabel('n'); ylabel('Modularity Qb');
savefig('modularity_Qb_vs_n_random')
% Nestedness
figure(6); 
errorbar(para,nest(:,inet),nest_std(:,inet))
xlabel('n'); ylabel('Nestedness');
savefig('nestedness_vs_n_random')
% Interia/exteria edge ratio
figure(7);
errorbar(para,ierat(:,inet),ierat_std(:,inet))
xlabel('n'); ylabel('Interia to exteria links ratio');
savefig('modularity_Qr_vs_n_random')