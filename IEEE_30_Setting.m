% 用来设置IEEE30节点下的参数
function [Setting,Para] = IEEE_30_Setting(Setting,xlsname,sheetname)

    Setting.NumP = 5;
    Setting_suffix = [num2str(Setting.NumP), 'P'];
    Setting.sample_method = 2; %0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
    Setting.cumulative_curve = 1;
    Setting.percent_strategic = 0.1;

    Setting.heuristic_spread_start = 0;
    Setting.heuristic_spread_end = 10;
    Setting.heuristic_spread_tolerance = 0.01;
    [Para,Num] = read_info([xlsname,'.xlsx'],sheetname, Setting);
    %% 考虑拓扑

    Setting.topology = 1;
    if Setting.topology == 1
        %导入拓扑的程序
        mpc = case30;
        GSDF_lb = makePTDF(mpc); %node和Branch之间的转换矩阵, Num.Branch *Num.Node.
        Para.GSDF = GSDF_lb(:,Para.node); %主体和Branch之间的转换矩阵, Num.Branch *Num.Participant.
        Para.GSDF(abs(Para.GSDF) < 1e-2) = 0;
        Para.BranchMax = mpc.branch(:,6);
        Para.BranchMax = Para.BranchMax/2;
    %     Para.BranchMax(find(Para.BranchMax == 0)) = 1000;
        Para.BranchMax = Para.BranchMax(sum(abs(Para.GSDF),2) > 0);
        Para.GSDF = Para.GSDF(sum(abs(Para.GSDF),2) > 0,:);
        Num.Branch = size(Para.BranchMax,1);
    end

    %%
    Para.SW_max = 0;%这个离散化了以后会有很多问题
    Setting.Qipp = 0;
    Setting.SWtolerance = 1e-5;
    Setting.Complementary_bigM = max(abs(min(Para.qmin)),abs(max(Para.qmax)));
    if Setting.topology == 1
        Setting.Complementary_bigM = max(Setting.Complementary_bigM, max(Para.BranchMax) * 2) * 1.1;
    end
    Setting.bigM = 100000;
    Setting.unified_spread = 1; % 是否设计统一的spread
    Setting.riskSW = 1;         % 是否设计带有风险项的SW计算
    Setting.direct_solve = 0;
    Setting.Pareto_Point = 21;
end
