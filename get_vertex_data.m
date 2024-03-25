function vertex = get_vertex_data(Result_vertex, Result_maxSW, Para, Setting, Num)
    % 整理vertex
    % vertex的四项分别是feasible,SW,BB,IR.IC,如果有第六项的话是maxSW
    % 展现各场景下的松弛量
    vertex.Relax = [[Result_vertex.RelaxSW]' [Result_vertex.RelaxIR]' ...
                [Result_vertex.RelaxBB]' [Result_vertex.RelaxIC]'];
    vertex.Relax(:,1) = vertex.Relax(:,1)/Para.SW_max;
    vertex.Relax(:,3) = vertex.Relax(:,3)/Para.SW_max;

    % 展现各场景下的Para
    vertex.para = [[Result_vertex.unified_spread]'*2  ...
                [Result_vertex.kIR]' [Result_vertex.kIC]'];

    % 计算QXPara
    for c = 1:Num.C %这个可以补充进来
        tempno = Para.company_set(c).no;
        Result_maxSW.QXPara(c) = Para.xpara(tempno) * Result_maxSW.QX(tempno);       
    end
    Result_maxSW.QXPara = Result_maxSW.QXPara(:);

    for i = 1:length(Result_vertex)
        for c = 1:Num.C
            tempno = Para.company_set(c).no;
            Result_vertex(i).QXPara(c) = Para.xpara(tempno) * Result_vertex(i).QX(tempno);       
        end
        Result_vertex(i).QXPara = Result_vertex(i).QXPara(:);
    end
    %计算LMPbalance
    Result_maxSW.LMPsettleC_RX = Result_maxSW.QXPara - Result_maxSW.LMPsettleC;
    Result_maxSW.LMPbalance = sum(Result_maxSW.LMPsettleC)-Result_maxSW.congestion_revenue;
    % 计算producer_welfare和consumer_welfare
    for i = 1:length(Result_vertex)
        Result_vertex(i).producer_welfare = sum(Result_vertex(i).RX(Para.Gset_company));
        Result_vertex(i).consumer_welfare = sum(Result_vertex(i).RX(Para.Dset_company));
    end
    Result_maxSW.producer_welfare = sum(Result_maxSW.LMPsettleC_RX(Para.Gset_company)); % 如果MaxSW不是用的基于LMP的结算，那其实对应的也不是RelaxIC
    Result_maxSW.consumer_welfare = sum(Result_maxSW.LMPsettleC_RX(Para.Dset_company));

    
    
    % 各场景下的welfare.  为什么maxSWde congestion revenue不一样？因为人家是真的spread =
    % 0,而我们这个spread实际是0.02.
    vertex.welfare = [[Result_vertex.qtotal]' ...
                [Result_vertex.welfare]' [Result_vertex.surplus]' [Result_vertex.baseSW]' ...
                [Result_vertex.congestion_revenue]'...
                [Result_vertex.producer_welfare]' [Result_vertex.consumer_welfare]'];
    
    % 这里用的都是LMP_settle，而不是用的是基于RX的Settle
    maxSW_welfare = [Result_maxSW.qtotal Result_maxSW.welfare ...
        Result_maxSW.LMPbalance Result_maxSW.baseSW Result_maxSW.congestion_revenue ...
        Result_maxSW.producer_welfare Result_maxSW.consumer_welfare];
    
    vertex.welfare = [vertex.welfare; maxSW_welfare];

    % 展示LMP和LMPI
    vertex.LMP_matrix = [Result_vertex.LMP]'; %这个LMP到底是LMP还是LMPI?
    vertex.LMP_matrix = [vertex.LMP_matrix;Result_maxSW.LMP']

    vertex.LMPI_matrix = [Result_vertex.LMPI]'; %这个LMP到底是LMP还是LMPI?
    vertex.LMPI_matrix = [vertex.LMPI_matrix;Result_maxSW.LMPI']

    vertex.RX_matrix = [Result_vertex.RX]';
    vertex.MX_matrix = [Result_vertex.MX]';
    vertex.Pl_matrix = [Result_vertex.Pl]';
    vertex.QX_matrix = [Result_vertex.QX]';

    sum(vertex.MX_matrix,2)
    sum(vertex.RX_matrix,2)
    sum(vertex.MX_matrix,2) - [Result_vertex.congestion_revenue]'
    sum(vertex.RX_matrix,2) + [Result_vertex.congestion_revenue]'

    vertex.QXPara_matrix = [Result_vertex.QXPara]'; % 这里算出来的根本就不是综合SW，而是某种场景下的SW。
    sum(vertex.QXPara_matrix,2)
    % 为什么最后这个问题平衡不了呢？感觉很奇怪.baseSW和SW不是一个东西
    
    % 这个又是用来干嘛的？
    vertex.MX_LMPrelaxSW = [Result_maxSW.LMPsettleC';Result_vertex(2).MX']; %原来的LMP结算，现在的最小成交价差里的结算
    vertex.RX_LMPrelaxSW = [Result_maxSW.LMPsettleC_RX';Result_vertex(2).RX']; %原来的LMP结算，现在的最小成交价差里的结算

    vertex.RX_matrix = [vertex.RX_matrix;Result_maxSW.LMPsettleC_RX']; % 考虑到与基础场景的对比
    vertex.MX_matrix = [vertex.MX_matrix;Result_maxSW.LMPsettleC']; % 考虑到与基础场景的对比
    vertex.MX_matrix = [vertex.MX_matrix;vertex.MX_LMPrelaxSW(1,:)]; % 考虑到与基础场景的对比
    % vertex.QX_matrix = [Result_vertex.QX];
    % vertex.Pl_matrix = [Result_vertex.Pl];
    % vertex.RX_matrix = [Result_vertex.Pl];
end