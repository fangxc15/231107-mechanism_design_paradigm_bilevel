% 此时的MX应该不苏子和
% MX = QX*X - \int Q(t) dt
% RX项
% MX项
% 100(2), 200(1), 300(1)
set(groot,'defaultfigurePosition',[400 400 660 450]);
set(0,'defaultfigurecolor','w'); %设置背景颜色为白色
set(groot,'defaultAxesFontWeight','bold');
set(groot,'defaultLegendFontSize',12);
set(groot,'defaultAxesFontSize',14);
% legend('boxoff');
X = [100,200,200,300,    300,400];
QX =[2,  2,     1,  1,     0,  0];
RX =[300,100,100,0,         0,0];
MX =[500,500,300,300,     0,  0];
MXLMP = [200,400,200,300,0,0];
RX =[150,50,50,0,         0,0];
MXrelax_half = [350,450,250,300,0,0];
%%
figure
[ax,p1,p2] = plotyy(X,QX, X,[MXLMP;MX],'plot','plot');
set(ax(1),'XColor','k','YColor','b'); %设置x轴为黑色，左边y（也就是y1）轴为蓝色
set(ax(2),'XColor','k','YColor','r'); %设置x轴为黑色，右边y（也就是y2）轴为红色

xlabel('申报成本,元/MWh');
ylabel(ax(1),'中标量,MW');
ylabel(ax(2),'发电收入,元');
set(ax,'Position',[0.16,0.16,0.65,0.7])
title('Q_A(X)和M_A(X)的关系');
set(p1,'linestyle','-','color','b');
set(p2(1),'linestyle','- -','color','g');
set(p2(2),'linestyle','- -','color','r');
h = legend('中标量函数Q(X)',"边际定价下的收入函数M'(X)",'激励相容的收入函数M(X)','Location','Best','NumColumns',1);
legend('boxoff')
set(ax(1),'ylim',[0,3]);
set(ax(2),'ylim',[0,600]);
grid on 
set(ax(1),'yTick',[0:0.5:3]);
set(ax(2),'yTick',[0:100:600]);
% exportgraphics(gca,'演示性Q_1(X)和M_1(X)的关系_1..jpg','Resolution',300)
print('-dpng','-r1000','演示性Q_1(X)和M_1(X)的关系_1.png');
% saveas(gcf,'演示性Q_1(X)和M_1(X)的关系_1.png','png','-r1000')
%%
figure
[ax,p1,p2] = plotyy(X,QX, X,[MXLMP;MX;MXrelax_half],'plot','plot');
set(ax(1),'XColor','k','YColor','b'); %设置x轴为黑色，左边y（也就是y1）轴为蓝色
set(ax(2),'XColor','k','YColor','r'); %设置x轴为黑色，右边y（也就是y2）轴为红色
set(ax,'Position',[0.16,0.16,0.7,0.7])

xlabel('申报成本,元/MWh');
ylabel(ax(1),'中标量,MW');
ylabel(ax(2),'发电收入,元');
title('Q_A(X)和M_A(X)的关系');
set(p1,'linestyle','-','color','b');
set(p2(1),'linestyle','- -','color','g');
set(p2(2),'linestyle','- -','color','r');
set(p2(3),'linestyle','- -','color','k');

% set(p2(1),'color',[0 0 0]);
legend('中标量函数Q(X)',"边际定价下的收入函数M'(X)",'激励相容的收入函数M(X)', ...
"0.5-激励相容的收入函数M''(X)")
legend('boxoff')

% set(ax(2),'legend',['福利函数R(X)','发电收入函数M(X)'])
set(ax(1),'ylim',[0,3]);
set(ax(2),'ylim',[0,600]);
grid on 
set(ax(1),'yTick',[0:0.5:3]);
set(ax(2),'yTick',[0:100:600]);
print('-dpng','-r1000','演示性Q_1(X)和M_1(X)的关系_2.png');

% saveas(1,'演示性Q_1(X)和M_1(X)的关系_2.jpg','resolution',1000)


% Npoint = 21;
% x = linspace(0,10,Npoint);
% y1 = sin(x);
% y2 = cos(x);
% figure
% [ax,p1,p2] = plotyy(x,y1,x,y2,'plot','plot');
% set(ax(1),'XColor','k','YColor','b'); %设置x轴为黑色，左边y（也就是y1）轴为蓝色
% set(ax(2),'XColor','k','YColor','r'); %设置x轴为黑色，右边y（也就是y2）轴为红色
% xlabel('X');
% ylabel(ax(1),'Y1'); % left y-axis
% ylabel(ax(2),'Y2'); % right y-axis
% title('The relationship of X, Y1 and X, Y2')
% %以下两行分别设置数据y1,y2的线型、颜色、填充点类型及颜色
% set(p1,'linestyle','-','marker','o','color','b','MarkerFaceColor','b');
% set(p2,'linestyle','- -','marker','^','color','r','MarkerFaceColor','r');
% %以下两行分别设置v数据y1,y2的间距大小
% set(ax(1),'yTick',[-1:0.2:1]);
% set(ax(2),'yTick',[-1:0.1:1]);

