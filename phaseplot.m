%% Plotting the phase diagrams
% Modularity
figure(1);
imagesc(modu); 
xlabel('pV'); ylabel('pH');
xticklabels = 0:0.1:1;
xticks = linspace(1, size(modu, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:0.1:1;
yticks = linspace(1, size(modu, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'YDir', 'normal')
colormap jet;
colorbar;
savefig('modularity_Qb')
% Nestedness
figure(2);
imagesc(nest); 
xlabel('pV'); ylabel('pH');
xticklabels = 0:0.1:1;
xticks = linspace(1, size(nest, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:0.1:1;
yticks = linspace(1, size(nest, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'YDir', 'normal')
colormap jet; 
colorbar;
savefig('nestedness')
% Interia/exteria edge ratio
figure(3);
imagesc(ierat);
xlabel('pV'); ylabel('pH');
xticklabels = 0:0.1:1;
xticks = linspace(1, size(ierat, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:0.1:1;
yticks = linspace(1, size(ierat, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'YDir', 'normal')
colormap jet;
colorbar;
savefig('modularity_Qr')
% Aspect ratio of infection network
figure(4);
imagesc(asp);
xlabel('pV'); ylabel('pH');
xticklabels = 0:0.1:1;
xticks = linspace(1, size(ierat, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
yticklabels = 0:0.1:1;
yticks = linspace(1, size(ierat, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels, 'YDir', 'normal')
colormap jet;
colorbar;
savefig('aspect_ratio')