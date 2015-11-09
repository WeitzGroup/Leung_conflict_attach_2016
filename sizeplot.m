%% Plotting the dependence of modularity/nestedness on network size
% Modularity
figure(1);
plot(para,modu)
xlabel('n'); ylabel('Modularity Qb');
savefig('modularity_Qb_vs_n')
% Nestedness
figure(2); 
plot(para,nest)
xlabel('n'); ylabel('Nestedness');
savefig('nestedness_vs_n')
% Interia/exteria edge ratio
figure(3);
plot(para,ierat)
xlabel('n'); ylabel('Interia to exteria links ratio');
savefig('modularity_Qr_vs_n')