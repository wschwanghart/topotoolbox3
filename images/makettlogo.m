function makettlogo


cols = [1 1 1;...
        0 0.4470 0.7410;...
        0.8500 0.3250 0.0980;...
        0.9290 0.6940 0.1250];


fn = 'Helvetica';

w = 0.5;
h = w;
va = 'middle';
ha = 'center';

figure()
ax = gca;
axis(ax,"square");
ylim(ax,[-h/2 -h/2+h ])
xlim(ax,[-w/2 -w/2+w ])

hrect = rectangle('Position',[-w/2 -h/2 w h],'Curvature',0.05,'LineWidth',15,...
    'FaceColor',cols(2,:),'EdgeColor',cols(2,:));

n = 6.5;
% ax = gca;
% ax.Position = [-w/2 -h/2 w h];


text(0,0,' T  ','FontSize',48*n,'FontWeight','bold','Color',cols(4,:),...
    'VerticalAlignment',va,'HorizontalAlignment',ha,...
    'FontUnits','inches','FontName',fn)
text(0,0,' T','FontSize',48*n,'FontWeight','bold','Color',cols(3,:),...
    'VerticalAlignment',va,'HorizontalAlignment',ha,...
    'FontUnits','inches','FontName',fn)
text(-0.02,-0.045,'  3','FontSize',44*n,'FontWeight','bold','Color',cols(1,:),...
    'VerticalAlignment',va,'HorizontalAlignment',ha,...
    'FontUnits','inches','FontName',fn)

set(gca,'XColor','none','YColor','none')
