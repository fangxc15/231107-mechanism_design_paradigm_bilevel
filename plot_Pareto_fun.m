function [] = plot_Pareto_fun(ifsave, Para, Result_ParetoBBSW,...
                    Result_ParetoBBIR, Result_ParetoBBIC,Result_ParetoSWIR,...
                    Result_ParetoSWIC, Result_ParetoIRIC, Result_ParetoSWIRBB, Result_primal_collect,...
                    vertex, Picture_folder, varargin)
    if length(varargin) >=1
        iftitle = cell2mat(varargin(1));
    else
        iftitle = 1;
    end
                
    Num.C = size(vertex.RX_matrix,2);
                
    set(groot,'defaultfigurePosition',[200 200 340 280]);
%     set(groot,'defaultLegendFontName','Times New Roman');
    set(groot,'defaultLegendFontSize',14);
    set(groot,'defaultAxesFontSize',14);
%     set(groot,'defaultFontSize',14);

    set(groot,'defaultAxesFontWeight','bold');
%     set(groot,'defaultAxesFontName','Times New Roman');
    set(groot,'defaultAxesFontName',['SimSun']);
    set(0,'defaultfigurecolor','w'); %设置背景颜色为白色
% %     ifsave = 0;             
%     version_suffix =  '_newversion2';            
%     Picture_root_folder = ['Picture', version_suffix];  
%     mkdir(Picture_root_folder);
%     Picture_folder = [Picture_root_folder,'/',save_name];
%     mkdir(Picture_folder);
    %% 画出BB,SW的帕累托前沿(1)
    
    figure(1)
    plotmatrix = sortrows([[Result_ParetoBBSW.welfare]' [Result_ParetoBBSW.surplus]' [Result_ParetoBBSW.kIR]'],1);
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
%     matrix_plot = sortrows([[Result_ParetoBBSW(solveindex).surplus]
%     solveindex = find([solution_ParetoBBSW.problem] == 0); % 要找到有解的,所以soluton问题要查看是怎么回事
%     [Result_ParetoBBSW(solveindex).welfare]]',1);
%     plot(matrix_plot(:,1),matrix_plot(:,2));
    grid on
    xlabel('社会福利 ($)');
    ylabel('收支盈余 ($)');
    if iftitle
        title('关于社会福利和盈余的帕累托曲线')
    end
    if ifsave
%         print('-dpng','-r1000',[Picture_folder,'/','BB_SW_Pareto.png']);
        saveas(1,[Picture_folder,'/','BB_SW_Pareto.jpg'])
    end
    %% 画出RelaxSW,RelaxBB的帕累托前沿(1)
    figure(2)
    plotmatrix = sortrows([[Result_ParetoBBSW.RelaxSW]' [Result_ParetoBBSW.RelaxBB]' ...
        [Result_ParetoBBSW.unified_spread]' [Result_ParetoBBSW.kIR]' [Result_ParetoBBSW.kIC]'],1);
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    
    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plotmatrix(:,2) = plotmatrix(:,2)/Para.SW_max;
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
    plotSWBB = plotmatrix;
    xlabel('社会福利的松弛量');
    ylabel('收支平衡的松弛量');
    grid on
    ylim([-0.5,0.5]);
    set(gca,'ytick',-0.5:0.25:0.5);
    if iftitle
        title('关于社会福利松弛项和收支平衡松弛项的帕累托曲线')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','Relax_BB_SW_Pareto.png']);
        saveas(2,[Picture_folder,'/','Relax_BB_SW_Pareto.jpg'])
    end
    %% RelaxSW, RelaxBB (这个不是前沿，而是离散点，来自Result_primal_collect)
    
    plotmatrix = sortrows([[Result_primal_collect.RelaxSW]' [Result_primal_collect.RelaxBB]' ...
        [Result_primal_collect.unified_spread]' [Result_primal_collect.kIR]' ...
        [Result_primal_collect.kIC]'],1);
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    
    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plotmatrix(:,2) = plotmatrix(:,2)/Para.SW_max;
%     plotmatrix(:,3) = plotmatrix(:,3) * 2;
    % scatter的点
    figure(3)
    scatter(plotmatrix(:,1), plotmatrix(:,2))
    
    grid on
    xlabel('社会福利的松弛量');
    ylabel('收支平衡的松弛量');
%     grid on
%     ylim([-0.5,0.5]);
%     set(gca,'ytick',-0.5:0.25:0.5);
    if iftitle
        title('关于社会福利松弛项和收支平衡松弛项的离散点')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','Relax_BB_SW_Pareto_else.png']);
        saveas(3,[Picture_folder,'/','Relax_BB_SW_Pareto_else.jpg'])
    end
    % 随着spread如何变化
    figure(4)
    plotmatrix = sortrows(plotmatrix,3);
    plot(plotmatrix(:,3), plotmatrix(:,1:2),'Linewidth',2)
    grid on
    xlabel('最小成交价差 ($/MWh)');
    ylabel('收支平衡/社会福利的松弛量');
%     grid on
%     ylim([-0.5,0.5]);
%     set(gca,'ytick',-0.5:0.25:0.5);
    if iftitle
        title('收支平衡/社会福利的松弛量随最小成交价差的变化')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','Relax_BB_SW_spreadfunc.png']);
        saveas(4,[Picture_folder,'/','Relax_BB_SW_spreadfunc.jpg'])
    end

    %% 画出BB,IR的帕累托前沿
    figure(5)
    plotmatrix = sortrows([[Result_ParetoBBIR.RelaxBB]' [Result_ParetoBBIR.RelaxIR]' ...
        [Result_ParetoBBIR.kIR]'],1);
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    %为什么会出现直线问题? 因为本身的出清结果是已经固定的。由于边际上的RX没有冲抵的区间，所以只能这样。
    %修改IR的定义方式后解决了该问题
    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
    plotBBIR = plotmatrix;
    grid on
    set(gca,'xtick',-0.5:0.25:0.5);
    xlabel('收支平衡松弛量');
    ylabel('个体理性松弛量');
    if iftitle
        title('关于收支平衡和个体理性松弛量的帕累托曲线')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','BB_IR_Pareto.png']);
        saveas(5,[Picture_folder,'/','BB_IR_Pareto.jpg'])
    end
    
    %% 画出BB,IC的帕累托前沿
    figure(6)
    plotmatrix = sortrows([[Result_ParetoBBIC.RelaxBB]' [Result_ParetoBBIC.RelaxIC]' ...
        [Result_ParetoBBIC.kIC]'],1);
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    %为什么会出现直线问题? 因为本身的出清结果是已经固定的。由于边际上的RX没有冲抵的区间，所以只能这样。
    %修改IR的定义方式后解决了该问题
    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plotmatrix(:,2) = plotmatrix(:,2)* min(Para.xpara);
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
    plotBBIC = plotmatrix;
    grid on
    set(gca,'xtick',-0.5:0.25:0.5);
    xlabel('收支平衡松弛量');
    ylabel('激励相容松弛量');
    if iftitle
        title('关于收支平衡和激励相容松弛量的帕累托曲线')
    end
        
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','BB_IC_Pareto.png']);
        saveas(6,[Picture_folder,'/','BB_IC_Pareto.jpg'])
    end
  %% 画出SW,IR的帕累托前沿
    figure(7)
    plotmatrix = [[Result_ParetoSWIR.RelaxSW]' [Result_ParetoSWIR.RelaxIR]' ...
        [Result_ParetoSWIR.kIR]'];

    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    %为什么会出现直线问题? 因为本身的出清结果是已经固定的。由于边际上的RX没有冲抵的区间，所以只能这样。
    %修改IR的定义方式后解决了该问题
    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
    plotSWIR = plotmatrix;
    grid on
%     set(gca,'xtick',-0.5:0.25:0.5);
    xlabel('社会福利松弛量');
    ylabel('个体理性松弛量');
    if iftitle
        title('关于社会福利和个体理性松弛量的帕累托曲线')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','SW_IR_Pareto.png']);
        saveas(7,[Picture_folder,'/','SW_IR_Pareto.jpg'])
    end
    %% 画出SW,IC的帕累托前沿
    figure(8)
    plotmatrix = [[Result_ParetoSWIC.RelaxSW]' [Result_ParetoSWIC.RelaxIC]'];

    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    %为什么会出现直线问题? 因为本身的出清结果是已经固定的。由于边际上的RX没有冲抵的区间，所以只能这样。
    %修改IR的定义方式后解决了该问题
    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plotmatrix(:,2) = plotmatrix(:,2)* min(Para.xpara);
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
    plotSWIC = plotmatrix;
    grid on
%     set(gca,'xtick',-0.5:0.25:0.5);
    xlabel('社会福利松弛量');
    ylabel('激励相容松弛量');
    if iftitle
        title('关于社会福利和激励相容松弛量的帕累托曲线')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','SW_IC_Pareto.png']);
        saveas(8,[Picture_folder,'/','SW_IC_Pareto.jpg'])
    end

    %% 画出IR,IC的帕累托前沿
    plotmatrix = [[Result_ParetoIRIC.RelaxIR]' [Result_ParetoIRIC.RelaxIC]' ...
                [Result_ParetoIRIC.kIR]' [Result_ParetoIRIC.kIC]'];
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
    % plotmatrix = plotmatrix(find([solution_ParetoIRIC.problem] == 0),:);
    plotmatrix = sortrows(plotmatrix,1);
    %在RelaxBB>0的时候，此时也能够Relax部分IC, 因为此时也能让k^IR有部分补偿
    %这里应该要把部分infeasible去除掉
    figure(9)
    plotmatrix(:,2) = plotmatrix(:,2)* min(Para.xpara);
    plot(plotmatrix(:,1), plotmatrix(:,2),'Linewidth',2)
    plotIRIC = plotmatrix;
    grid on
%     set(gca,'xtick',-0.5:0.25:0.5);
    xlabel('个体理性松弛量');
    ylabel('激励相容松弛量');
    if iftitle
        title('关于个体理性和激励相容松弛量的帕累托曲线')
    end
    if ifsave
        print('-dpng','-r1000',[Picture_folder,'/','IR_IC_Pareto.png']);
        saveas(9,[Picture_folder,'/','IR_IC_Pareto.jpg'])
    end
%% 画出六联图
    set(groot,'defaultLegendFontSize',14);
    set(groot,'defaultAxesFontSize',15);
    figure(21)
    set(gcf,'Position',[0 0 1200 1000]);
    
    subplot(3,2,1)
    h1 = plot(plotSWBB(:,2), plotSWBB(:,1),'LineWidth',2);
    set(gca,'YTick',0:0.03:0.12)
    set(gca,'XTick',-0.5:0.25:0.5)
    set(gca,'looseInset',[0 0 0 0])
    set(gca, 'OuterPosition',[0.01 0.677 0.48 0.313]);
    h1_axis = gca;
    xlabel('收支平衡松弛量');
    ylabel('社会福利松弛量');
    ylim([0, 0.12])
    xlim([-0.5,0.5])
    grid on
    
    subplot(3,2,3)
    h3 = plot(plotBBIC(:,1), plotBBIC(:,2),'LineWidth',2); 
    set(gca,'XTick',-0.5:0.25:0.5)
    set(gca,'YTick',0:0.15:0.6)
    set(gca,'looseInset',[0 0 0 0])
    set(gca, 'OuterPosition',[0.01,0.343,0.48,0.313]);
%     set(gca,'xtick',-0.5:0.25:0.5);
    h3_axis = gca;
    xlim([-0.5, 0.5])
    ylim([0,0.6])
    xlabel('收支平衡松弛量');
    ylabel('激励相容松弛量');
    grid on

    
    subplot(3,2,5)
    h5 = plot(plotBBIR(:,1), plotBBIR(:,2),'LineWidth',2);
    set(gca,'XTick',-0.5:0.25:0.5)
    set(gca,'YTick',0:0.05:0.2)
    set(gca,'looseInset',[0 0 0 0])
    set(gca, 'OuterPosition',[0.01,0.010,0.48,0.313]);
    h5_axis = gca;
    xlim([-0.5, 0.5])
    ylim([0,0.2])
    xlabel('收支平衡松弛量');
    ylabel('个体理性松弛量');
    grid on
    
    
    subplot(3,2,2)
    h2 = plot(plotSWIR(:,2), plotSWIR(:,1),'LineWidth',2);
    set(gca,'XTick',0:0.02:0.08)
    set(gca,'YTick',0:0.01:0.04)
    set(gca,'looseInset',[0 0 0 0])
    set(gca, 'OuterPosition',[0.51,0.677,0.48,0.313]);
    h2_axis = gca;
    xlim([0, 0.08])
    ylim([0,0.04])
%     set(gca,'xtick',-0.5:0.25:0.5);
    xlabel('个体理性松弛量');
    ylabel('社会福利松弛量');
    grid on
    
    
    subplot(3,2,4)
    h4 = plot(plotIRIC(:,1), plotIRIC(:,2),'LineWidth',2);
    set(gca,'XTick',0:0.02:0.08)
    set(gca,'YTick',0:0.05:0.2)
    set(gca,'looseInset',[0 0 0 0])
    set(gca, 'OuterPosition',[0.51,0.343,0.48,0.313]);
    h4_axis = gca;
    xlim([0,0.08])
    ylim([0,0.2])
    xlabel('个体理性松弛量');
    ylabel('激励相容松弛量');
    grid on

    
    
    subplot(3,2,6)
    h6 = plot(plotSWIC(:,1), plotSWIC(:,2),'LineWidth',2);
    set(gca,'XTick',0:0.01:0.04)
    set(gca,'YTick',0:0.05:0.2)
    set(gca,'looseInset',[0 0 0 0])
    set(gca, 'OuterPosition',[0.51,0.010,0.48,0.313]);
    h6_axis = gca;
    xlim([0,0.04])
    ylim([0,0.2])
    xlabel('社会福利松弛量');
    ylabel('激励相容松弛量');
    grid on

%     set(gca,'xtick',-0.5:0.25:0.5);


    

    if iftitle
        title('机制松弛项间的两两权衡')
    end

%     Expand_axis_fill_figure(h1_axis);
%     Expand_axis_fill_figure(h2_axis);
%     Expand_axis_fill_figure(h3_axis);
%     Expand_axis_fill_figure(h4_axis);
%     Expand_axis_fill_figure(h5_axis);
%     Expand_axis_fill_figure(h6_axis);
    
    
%     set(gca,'xtick',-0.5:0.25:0.5);
	if ifsave
        print('-dpng','-r600',[Picture_folder,'/','Pareto_collect.png']);
        saveas(21,[Picture_folder,'/','Pareto_collect.jpg'])
    end
%% 画出SWIRBB的Pareto Frontier
    set(groot,'defaultAxesFontSize',11);
    set(groot,'defaultLegendFontSize',11);
    plotmatrix = [[Result_ParetoSWIRBB.RelaxSW]' [Result_ParetoSWIRBB.RelaxIR]' ...
                [Result_ParetoSWIRBB.RelaxBB]' [Result_ParetoSWIRBB.unified_spread]'...
                [Result_ParetoSWIRBB.kIR]' [Result_ParetoSWIRBB.kIC]'];
    % plotmatrix = plotmatrix(find([solution_ParetoSWIRBB.problem] == 0),:);
    plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0 | ...
        plotmatrix(:,3) ~= 0 ,:);

    plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
    plotmatrix(:,3) = plotmatrix(:,3)/Para.SW_max;

    set(groot,'defaultfigurePosition',[300 300 460 390]);
    set(0,'defaultfigurecolor','w'); %设置背景颜色为白色
    set(groot,'defaultAxesFontWeight','bold');

    [X,Y] = meshgrid(linspace(min(plotmatrix(:,1)),max(plotmatrix(:,1)),20), ...
            linspace(min(plotmatrix(:,2)),max(plotmatrix(:,2)),20));
    
    Z = griddata(plotmatrix(:,1),plotmatrix(:,2),plotmatrix(:,3),X,Y);
    figure(10)
    s = surf(X,Y,Z);
    s.FaceAlpha = 0.9;
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    set(gca,'XDir','reverse')
    set(gca,'YDir','reverse')
    xlabel('社会福利松弛量');
    % ylim([0,0.08]);
    ylabel('个体理性松弛量');
    zlabel('收支平衡松弛量');
    if iftitle
        title('社会福利-个体理性-收支平衡松弛的帕累托前沿')
    end
    if ifsave
        print('-dpng','-r600',[Picture_folder,'/','SWIRBB_Pareto.png']);
        saveas(10,[Picture_folder,'/','SWIRBB_Pareto.jpg'])
    end
    %%
    set(groot,'defaultAxesFontSize',15);
    set(groot,'defaultLegendFontSize',15);
    colors = [0.8,0.56,0.4;%红色,RelaxSW
              0.22,0.80,0.57;%绿色,RelaxBB
              1,0.5,0; %橙色,RelaxIR
              0.5,0.7,0.9; %蓝色,RelaxIC
              0.5,0.7,0.9; %蓝色,LMP
              ];
    
    
%     [0.93 0.92 0.27; % 原来的社会福利最大化情况
%           0.16 0.87 0.58;     %
%           0.5 0.7 0.9];
% 
    
    %% 画图比较出清量
    set(groot,'defaultfigurePosition',[200 200 710 420]);

    figure(11)
    b = bar(vertex.QX_matrix([5 2],:)');
    xlabel('机组/用户编号')
    ylabel('出清电量 (MWh)')
    if iftitle
         title('松弛社会福利最大造成的出清量降低')
    end
    set(gca,'Ticklength',[0,0])

    legend('社会福利最大下的出清','松弛社会福利下的出清','Location','SouthEast');
    legend('boxoff')
    set(b(1),'FaceColor',colors(4,:))     % 颜色调节
    set(b(2),'FaceColor',colors(1,:))     % 颜色调节
    set(gca,'YGrid','on')

    if ifsave
        print('-dpng','-r600',[Picture_folder,'/','QX_compare.png']);
        saveas(11,[Picture_folder,'/','QX_compare.jpg'])
    end

    %% 画图比较潮流量
    set(groot,'defaultfigurePosition',[200 200 1010 420]);

    figure(12)
    b = bar(vertex.Pl_matrix([5 2],:)');
    xlabel('线路编号')
    ylabel('潮流量 (MWh)')
    if iftitle
        title('松弛社会福利最大造成的线路潮流变化')
    end
    legend('社会福利最大下的出清','松弛社会福利下的出清','Location','NorthEast');
    legend('boxoff')
    set(b(1),'FaceColor',colors(4,:)) % 颜色调节
    set(b(2),'FaceColor',colors(1,:))     % 颜色调节
    set(gca,'YGrid','on')
    set(gca,'Ticklength',[0,0])

    if ifsave
        print('-dpng','-r600',[Picture_folder,'/','Pl_compare.png']);
        saveas(12,[Picture_folder,'/','Pl_compare.jpg'])
    end

% set(gca,'yColor',[0.25,0.41,0.88],'FontWeight','bold'); %纵坐标轴颜色设置

    %% 要比较所有的场景的welfare
    set(groot,'defaultfigurePosition',[200 200 910 420]);

    figure(31)
    b = bar(vertex.RX_matrix([2:4 6],:)');
    Num.C = size(vertex.RX_matrix,2);
    for i = 1:Num.C-1
        line([i+0.5,i+0.5],[-200,1000],'color', 'k','LineStyle','--')
    end
    xlim([0.5,10.5])
    xlabel('主体编号')
    ylabel('主体福利 ($)')
    if iftitle 
        title('各极点场景处的福利函数对比')
    end

    legend('松弛社会福利最大','松弛收支平衡(VCG)','松弛个体理性','松弛激励相容(LMP)','Location','NorthEast');
    legend('boxoff')
    set(b(1),'FaceColor',colors(1,:))     % 颜色调节
    set(b(2),'FaceColor',colors(2,:)) % 颜色调节
    set(b(3),'FaceColor',colors(3,:))     % 颜色调节
    set(b(4),'FaceColor',colors(5,:))     % 颜色调节
    set(gca,'Ticklength',[0,0])


    set(gca,'YGrid','on')

    if ifsave
        print('-dpng','-r600',[Picture_folder,'/','RX_compare_all.png']);
        saveas(31,[Picture_folder,'/','RX_compare_all.jpg'])
    end

%     %% 画图比较RX (VCG,RelaxSW)
%     set(groot,'defaultfigurePosition',[200 200 810 420]);
% 
%     figure(13)
%     b = bar(vertex.RX_matrix([3 2],:)');
%     xlabel('主体编号')
%     ylabel('效用量')
%     title('松弛社会福利最大造成的主体效用变化')
%     legend('社会福利最大出清(VCG)','松弛社会福利下的出清','Location','NorthEast');
%     legend('boxoff')
%     set(b(1),'FaceColor',colors(3,:)) % 颜色调节
%     set(b(2),'FaceColor',colors(2,:))     % 颜色调节
%     set(gca,'YGrid','on')
%     if ifsave
%         print('-dpng','-r600',[Picture_folder,'/','RX_compare.png']);
%         saveas(13,[Picture_folder,'/','RX_compare.jpg'])
%     end
%     %% 要画出原本的LMP机制和现在的引入min_spread后的福利函数区别
%     
%     set(groot,'defaultfigurePosition',[200 200 610 420]);
%     set(groot,'defaultLegendFontSize',16);
% %     set(groot,'defaultAxesFontSize',11);
% %     set(groot,'defaultFontSize',14);
%     figure(17)
%     b = bar(vertex.RX_LMPrelaxSW');
%     xlabel('主体编号')
%     ylabel('主体的支付总额')
%     title('松弛社会福利后的福利变化')
%     legend('原情况节点电价下的福利','松弛社会福利场景下的福利','Location','NorthEast');
%     legend('boxoff')
%     set(gca,'YGrid','on')
%     set(b(1),'FaceColor',colors(1,:)) % 颜色调节
%     set(b(2),'FaceColor',colors(2,:))     % 颜色调节
%     if ifsave
%         print('-dpng','-r600',[Picture_folder,'/','RX_compare_LMPnew_minspread.png']);
%         saveas(17,[Picture_folder,'/','RX_compare_LMPnew_minspread.jpg'])
%     end
%     %% 要画出RelaxIC, RelaxIR, RelaxBB 三种情况下不同的welfare分布
%     set(groot,'defaultfigurePosition',[200 200 810 420]);
% 
%     figure(14)
%     b = bar(vertex.RX_matrix([6 3:4],:)');
%     xlabel('主体编号')
%     ylabel('主体福利')
%     title('松弛BB,IR,IC场景下的主体福利对比')
%     legend('松弛激励相容(LMP)的社会福利','松弛收支平衡(VCG)下的主体福利','松弛个体理性下的社会福利','Location','NorthEast');
%     legend('boxoff')
%     set(b(1),'FaceColor',colors(1,:))     % 颜色调节
%     set(b(2),'FaceColor',colors(3,:)) % 颜色调节
%     set(b(3),'FaceColor',colors(4,:))     % 颜色调节
%     set(gca,'YGrid','on')
% 
%     if ifsave
%         print('-dpng','-r600',[Picture_folder,'/','RX_compare_RelaxBBIRIC.png']);
%         saveas(14,[Picture_folder,'/','RX_compare_RelaxBBIRIC.jpg'])
%     end

    %% 要比较所有的结算函数区别
    set(groot,'defaultfigurePosition',[200 200 910 420]);
    set(groot,'defaultLegendFontSize',16);
    
%     vertex.MX_matrix = [vertex.MX_matrix;vertex.MX_LMPrelaxSW(1,:)]; % 考虑到与基础场景的对比

    figure(32)
    b = bar([vertex.MX_matrix([2:4],:);vertex.MX_LMPrelaxSW(1,:)]');
    for i = 1:Num.C-1
        line([i+0.5,i+0.5],[-4000,3000],'color', 'k','LineStyle','--')
    end
    xlim([0.5,10.5])
    xlabel('主体编号')
    ylabel('主体结算值 ($)')
    if iftitle 
        title('各极点场景处的结算函数对比')
    end
    legend('松弛社会福利最大','松弛收支平衡(VCG)','松弛个体理性','松弛激励相容(LMP)','Location','NorthWest');
    legend('boxoff')
    set(b(1),'FaceColor',colors(1,:))     % 颜色调节
    set(b(2),'FaceColor',colors(2,:)) % 颜色调节
    set(b(3),'FaceColor',colors(3,:))     % 颜色调节
    set(b(4),'FaceColor',colors(5,:))     % 颜色调节
    set(gca,'Ticklength',[0,0])

    set(gca,'YGrid','on')

    if ifsave
        print('-dpng','-r600',[Picture_folder,'/','MX_compare_all.png']);
        saveas(32,[Picture_folder,'/','MX_compare_all.jpg'])
    end

    % 要画出原本的LMP机制和现在的引入min_spread后的结算函数区别
    
%     set(groot,'defaultAxesFontSize',11);
%     set(groot,'defaultFontSize',14);
%     figure(16)
%     b = bar(vertex.MX_LMPrelaxSW');
%     xlabel('主体编号')
%     ylabel('主体的支付总额')
%     title('松弛社会福利后的支付总额变化')
%     legend('原情况节点电价下的支付','松弛社会福利场景下的支付','Location','NorthWest');
%     legend('boxoff')
%     set(gca,'YGrid','on')
%     set(b(1),'FaceColor',colors(1,:)) % 颜色调节
%     set(b(2),'FaceColor',colors(2,:))     % 颜色调节
%     if ifsave
%         print('-dpng','-r600',[Picture_folder,'/','MX_compare_LMPnew_minspread.png']);
%         saveas(16,[Picture_folder,'/','MX_compare_LMPnew_minspread.jpg'])
%     end
    
    %% 要画出普通的LMP和松弛社会福利后的LMP变化(这个LMP变化有意义吗？感觉不是很有意义啊）
    set(groot,'defaultfigurePosition',[200 200 810 420]);
    figure(15)
    b = plot(vertex.LMPI_matrix([6,2],:)','LineWidth',2);
    b(1).Color = colors(1,:);
    b(2).Color = colors(2,:);
    xlabel('机组编号')
    ylabel('机组所在节点的节点电价 ($/MWh)')
    if iftitle
        title('引入最小成交价差后的节点电价变化')
    end
    legend('原情况下的节点电价','引入最小成交价差后的节点电价','Location','NorthWest');
    legend('boxoff')
    set(gca,'YGrid','on')
    set(gca,'Ticklength',[0,0])

%     set(b(1),'FaceColor',[0,0.78,0.55]) % 颜色调节
%     set(b(2),'FaceColor',[1,0.5,0])     % 颜色调节
    if ifsave
        print('-dpng','-r600',[Picture_folder,'/','LMP_compare_minspread.png']);
        saveas(15,[Picture_folder,'/','LMP_compare_minspread.jpg'])
    end

end
