%% Plotting the network
% Original matrix
figure(1); 
%set(gcf,'Position',[0 72 1751 922]);
bp.plotter.font_size = 2.0;
bp.plotter.PlotMatrix();
savefig('matrix')
% Nested matrix
figure(2); 
%set(gcf,'Position',[0+50 72 932 922]);
bp.plotter.use_isocline = true;
bp.plotter.isocline_color = 'red';
bp.plotter.PlotNestedMatrix(); 
savefig('nest_matrix')
% Plot modular network
figure(3);
bp.plotter.use_isocline = true;
bp.plotter.PlotModularMatrix();
title(['$Q = $',num2str(bp.community.Qb),' $c = $', num2str(bp.community.N)],...
'interpreter','latex','fontsize',23);
savefig('modu_matrix')
% Plotting in graph format
figure(4);
bp.plotter.PlotGraph();
savefig('graph')