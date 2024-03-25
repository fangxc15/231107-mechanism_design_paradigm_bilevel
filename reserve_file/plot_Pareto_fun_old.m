function [] = plot_Pareto_fun(ifsave, solution_ParetoBBSW,Result_ParetoBBSW,Result_ParetoBBIC,Result_ParetoIRSW,save_folder )
    mkdir(['Picture/',save_folder]);
    %% 画出BB,SW的帕累托前沿(1)
    solveindex = find([solution_ParetoBBSW.problem] == 0); % 要找到有解的,所以soluton问题要查看是怎么回事
    figure(1)
    matrix_plot = sortrows([[Result_ParetoBBSW(solveindex).surplus]
    [Result_ParetoBBSW(solveindex).welfare]]',1);
    plot(matrix_plot(:,1),matrix_plot(:,2));
    xlabel('盈余');
    ylabel('社会福利');
    title('关于社会福利和盈余的帕累托曲线')
    if ifsave
        saveas(1,['Picture/',save_folder,'/','BB_SW_Pareto.jpg'])
    end
    % 一般来讲，社会福利最大的点就对应着高低匹配的出清规则。此时会有收支不平衡，取决于每个主体创造的增量社会福利之和与社会总福利谁大
    % 如果前者大于后者,那么为了同时满足社会福利最大/激励相容
    % 损失一部分社会福利，可以换取收支平衡。进一步损失，可以获取收支盈余。但收支盈余有个极限，因为随着社会福利的减少，收支盈余量本来的潜力也在减少。

    % welfare的崩溃也有个极限.如果welfare再降低促进不了盈余量的话，它不会在帕累托前沿上.
    % 最后一点welfare的增加,需要付出极大的收支不平衡量.

    % 将市场从q(x)过渡到规则. 但这里自由度太高了，未必能获得一个合适的规则。是否可以去逼近这个规则? 建立一个映射. 这是一个好思路
%     %% 画出随着spread增加时候点的变化. 此时一直是保持个体理性和激励相容.通过修改q(x)以及对应的
%     % 这是规则化的点，但不一定在帕累托前沿上,帕累托前沿上的点也不能规则化。
%     matrix_plot_spread = [[result_collect_spread.choose_spread];[result_collect_spread.welfare];[result_collect_spread.balance]];
%     figure(2)
%     plot(matrix_plot_spread(1,:),matrix_plot_spread([2:3],:));
%     xlabel('选择的最小成交价差');
%     ylabel('社会福利/盈余量');
%     legend('社会福利','盈余量');
%     title('不同最小成交价差下的机制情况');
%     if ifsave
%         saveas(2,['Picture/',save_folder,'/','spread_BBSW.jpg'])
%     end
%     %可以看到,如果让最小成交价差为负,既会有损社会福利又会有损盈余量.
%     %由于r(x)是q(x)的积分,所以成交量越多需要付出分配给主体的效用必然越大,如果没有切实增大社会总福利,那么收支不平衡量必然越来越大.
%     %在一定的程度下，如果增大社会福利的速度快于付出效用的速度，那么收支盈余量将增大.但如果增大社会福利的速度慢于需要分配的主体效用,那么盈余量会减小.
%     %上图中,当价差从1逐步下降到0.25的时候,社会福利增加,盈余量也增加.但从0.25下降到0的时候,社会福利增加，盈余量减少，这说明社会福利增加的速度不够快了.
%     %所以帕累托前沿,应该是在[0,0.25]这一块.
%     %% 画出该情况下的帕累托前沿
%     figure(3)
%     plot(matrix_plot_spread(3,:),matrix_plot_spread(2,:));
%     hold on
%     plot(matrix_plot(:,1),matrix_plot(:,2))
%     xlabel('盈余量');
%     ylabel('社会福利');
%    legend('价差方式','全空间方式');
%     title('不同最小成交价差下的帕累托曲线');
%     
%     if ifsave
%         saveas(3,['Picture/',save_folder,'/','BB_SW_Pareto_2.jpg'])
%     end

    %% 画出BB,IC的帕累托前沿
    solveindex = 2:length(Result_ParetoBBIC);
    figure(4)
    matrix_plot = sortrows([[Result_ParetoBBIC(solveindex).RelaxBB]
    [Result_ParetoBBIC(solveindex).RelaxIC]]',1);
    plot(matrix_plot(:,1),matrix_plot(:,2));
    xlabel('收支平衡松弛量');
    ylabel('激励相容松弛量');
    title('关于收支平衡和激励相容松弛量的帕累托曲线')
    if ifsave
        saveas(4,['Picture/',save_folder,'/','BB_IC_Pareto.jpg'])
    end
  
    %% 画出IR,SW的帕累托前沿
    solveindex = 2:length(Result_ParetoIRSW);
    figure(5)
    matrix_plot = sortrows([[Result_ParetoIRSW(solveindex).RelaxIR]
    [Result_ParetoIRSW(solveindex).RelaxSW]]',1);
    plot(matrix_plot(:,1),matrix_plot(:,2));

    xlabel('个体理性松弛量');
    ylabel('社会福利松弛量');
    title('关于个体理性和社会福利松弛量的帕累托曲线')
    if ifsave
        saveas(5,['Picture/',save_folder,'/','IR_SW_Pareto.jpg'])
    end
end

