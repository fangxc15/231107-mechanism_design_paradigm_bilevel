clear 
script_type = 'xls';
xlsname = 'Clearing_IEEE30';
sheetname = 'Version1'; %这个地方开始出现不平衡的原因是把抽样方式改了。采用了间隔掺入式的抽样
folder_name = 'Result';

Setting.NumP = 5;
Setting_suffix = [num2str(Setting.NumP), 'P'];
Setting.sample_method = 2; %0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
Setting.cumulative_curve = 1;
Setting.percent_strategic = 0.1;


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

%% Heursitic find feasibile (find max SW, s.t. BB= 0) Optimization primal 
Result_maxSW = optimization_primal(Num,Para,Setting, zeros(1,Num.I));

Setting.percent_strategic = 1;
Result_LMP = optimization_primal(Num,Para,Setting, zeros(1,Num.I));
Setting.percent_strategic = 0.1;



Para.Dset_company = find(Para.sign == 1);
Para.Gset_company = find(Para.sign == -1);
sum(Result_LMP.RX(Para.Gset_company))
sum(Result_LMP.RX(Para.Dset_company))

Para.SW_max = Result_maxSW.SW;
Result_primal_BBbalance = heuristic_find_BBbalance(0,10,0.01,Num,Para,Setting);

% max(primal_SW(primal_BB >= 0))
Setting.Result_primal_SWmax = Result_primal_BBbalance(1); %社会福利最大的情况
Setting.bigM = Setting.Result_primal_SWmax.SW * 2;

primal_BB = [Result_primal_BBbalance.BB]';
primal_SW = [Result_primal_BBbalance.SW]';
primal_matrix = [primal_BB primal_SW (1:length(primal_BB))'];
primal_balance_matrix = sortrows(primal_matrix(primal_matrix(:,1) > 0,:),2,'descend');
BBbalanceno = primal_balance_matrix(1,3);
if length(primal_balance_matrix) >=1
    Setting.Result_primal_BBbalance = Result_primal_BBbalance(BBbalanceno); %保证BB的同时最大化社会福利的情况
end

% BBmaxno = find([Result_primal_BBbalance.BB] == max([Result_primal_BBbalance.BB]));
% Setting.Result_BBmax = Result_primal_BBbalance(BBmaxno);

%% Heuristic Find MaxBB
Result_primal_BBmax = heuristic_find_BBmax(0,10,0.01,Num,Para,Setting);

BBmaxno = find(([Result_primal_BBmax.BB] == max([Result_primal_BBmax.BB])));
Setting.Result_primal_BBmax = Result_primal_BBmax(BBmaxno(1));


%% Basic
tic;
[Var_basic, Cons_basic] = optimization_construct(Num,Para,Setting);
time.construct = toc;

%% Vertex
tic;
Setting.hard_assign = [0];

Para.BBIR_ratio = 1/sum(Para.Cstrategic); %Result_vertex(5).RelaxIR/Result_vertex(4).RelaxBB
for c = 1:Num.C
    tempno = Para.company_set(c).no;
    if Para.Cstrategic(c) == 1
        temp_len(c) = sum(sum(Para.Interval_len(tempno,:) .^2,1).^0.5);
    else
        temp_len(c) = 0;
    end
end         
Para.BBIC_ratio = 1/sum(temp_len); %Result_vertex(6).RelaxIC/Result_vertex(4).RelaxBB

[Result_vertex,solution_vertex] = solve_vertex(Setting, Para, Num, Var_basic, Cons_basic);
time.vertex = toc;

%% maxrelax
% maxrelax求解BB要比后续的Pareto前沿求解方便的多。因为此时是完全松弛掉了BB的约束，不会出现很麻烦的结果。
% 我之后的都是要满足BB的约束同时最大化SW（之前的满足BB平衡来最大化SW还是因为有启发式的解的缘故）
tic;
% Setting.Result_BBbalance = Result_vertex(3);
[Result_maxrelax,solution_maxrelax] = solve_maxrelax(Setting, Para, Num, Var_basic, Cons_basic);
time.maxrelax = toc;

for i = 1:length(Result_maxrelax)
    Result_maxrelax(i).RelaxSW = Result_vertex(4).SW - Result_maxrelax(i).SW;
end
%% More Heuristic Solutions
Result_primal_collect = Result_primal_BBbalance;
j = length(Result_primal_collect);
for i = 1:length(Result_primal_BBmax)
    j = j+1;
    Result_primal_collect(j) = Result_primal_BBmax(i);
end
temp_list = linspace(0,Result_maxrelax(1).unified_spread,30);
for i = 1:30
    j = j+1;
    Result_primal_collect(j) = optimization_primal(Num,Para,Setting, temp_list(i) * ones(1,Num.I));
end

Setting.primal_set = Result_primal_collect;
%%
% %% SW_BB
% % 这里是先给定一个SW的需求(relax SW的上限从小到大,SW的下限从大到小), 然后尝试最大化BB.
% % 这时候如果要满足SW,搜索方式应该很容易, 逐步增大spread即可。
% % 此时要最大化BB（或者最小化BBrelax),其实根据spread算出来的Q(X),配合S^{IC}, S^{IR}就可以得到MX, 进而求得BB
% % 之所以求解困难，是因为我由spread推导得到底层的KKT条件的时候，有很多整数。这些整数用来决策QX搭界或者mumax的情况，底下场景太多，也很难解。
% % 这里本质是整数太多，所以难以求解
% % 后面的根据S^IR, S^IC, Q(X)去直接决策M(X)，这个并不难
% % 核心点就在于如何去得到Q(X)和要求的SW。如果说SW要求的是最大的。
% % 
% % 基本上gap不降低，比BB_SW的降低速度还慢。本质上还是整数过多。本身决策Q(X)，但要求满足k的那一堆KKT，也很难。
% % 无论是SW，还是BB，其实本质都是要决策Q(X), 而Q(X)的决策都依赖互补松弛条件。
% % 
% % 这里的gap更大，问题更大。关键是，SW-BB法，从SW比较大的时候开始逐步减小，一开始SW的变化会有很大的BB变化。
% % 如果说我用BB-SW法，我一开始的BB比较大，然后逐步减小，当我BB比较大的时候，一开始BB的变化会导致很大的SW变化。
% % 
% % 所以理论上来讲我有两种改进的方式。
% % 一种是我在BB快到极限的时候，这个时候用SW等间距法，SW-BB去求，在SW快到极限的时候，用BB等间距法去求。
% % 另一种改进的方式是，赋予更好的启发式解法。
% 
% 
% Setting.Pareto_XVar = 'SW';
% Setting.Pareto_YVar = 'BB';
% Setting.set_IR = 0;
% Setting.set_IC = 0;
% % Result_vertex(4).welfare 代表最大的可供松弛的SW.
% 
% Setting.values = linspace(0,Result_vertex(4).SW - Result_maxrelax(1).SW,Setting.Pareto_Point);
% % 第二项表示最大的SW，减去BB最大情况下的SW
% 
% delete_point_num = round(Setting.Pareto_Point * 0.25) ;
% Setting.values = Setting.values(delete_point_num+1: Setting.Pareto_Point);
% % Setting.Result_primal_Paretostart = Result_vertex(4); % 一开始是社会福利最大的点?
% Setting.primal_set = Result_primal_collect; %Result_primal_BBmax(2:length(Result_primal_BBmax));
% 
% 
% % RelaxSW 一开始完全不松弛，后面松弛的越来越多, SW从大到小
% % linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);
% % RelaxBB, 一开始BB最大，后来BB越来越小.
% % 一开始的这几个点都需要舍弃的, 也就是RelaxSW太小的这几个要舍弃
% 
% % 
% % 一开始的求解肯定会很困难。我是不是要把它拼到一起比较好
% % 我写出了一个BB/SW的关系。应该是RelaxSW最小(SW最大)的这几个点舍弃掉，然后从一个RelaxSW比较小的点开始看起
% % RelaxSW比较的点，会对应一个还可以的BB。然后从这个BB，再向BB为负值的去看，得到剩余的一些。
% 
% 
% tic;
% [Result_ParetoSWBB,solution_ParetoSWBB] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
% time.ParetoSWBB = toc;
% Setting = rmfield(Setting, 'primal_set');
% Setting = rmfield(Setting, 'Result_primal_Paretostart');
% 
% 
% % BB_SW_temp (拼接上去的东西)
% Setting.Pareto_XVar = 'BB';
% Setting.Pareto_YVar = 'SW';
% Setting.set_IR = 0;
% Setting.set_IC = 0;
% 
% Setting.values = linspace(Result_ParetoSWBB(2).RelaxBB,...
%                         Result_vertex(4).RelaxBB, round(Setting.Pareto_Point * 0.75));
% Setting.primal_set = Result_primal_collect;
% 
% % 但按这样拼接有点蠢，应该反过来拼接。
%% BB_SW
% 最后的解不一定合理。比如这里，我把BB的阈值设计的很低了之后，我才发现有一个SW值在BB挺高的时候就能满足。那它显然是前几个Pareto
% Frontier上的可行解，但根本没被识别出来。
% 因此有的时候解就没搜索出来，感觉很尴尬
% 所以我的求解方式也很有讲究，怎么求解？现在是我的BB从大到小，前一个的解是后一个的可行解
% 如果我采用的方式是BB从小到大，那我前一个的解不一定是可行解，这样更难找。




% 说明是有gap, 没有求得最优解。
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'SW';
Setting.set_IR = 0;
Setting.set_IC = 0;
% Result_maxrelax(1).RelaxBB; %代表最大的可能盈余
% Result_vertex(4).RelaxBB;   %代表要求保证其他项时的盈余relax项
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);

delete_point_num = round(Setting.Pareto_Point * 0.25) ;
Setting.values = Setting.values(delete_point_num+1: Setting.Pareto_Point);



% Setting.Result_primal_Paretostart = Result_maxrelax(1);

tic;
[Result_ParetoBBSW,solution_ParetoBBSW] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBSW = toc;
if isfield(Setting,'primal_set')
    Setting = rmfield(Setting, 'primal_set');
end
if isfield(Setting,'Result_primal_Paretostart')
    Setting = rmfield(Setting, 'Result_primal_Paretostart');
end



Setting.Pareto_XVar = 'SW';
Setting.Pareto_YVar = 'BB';
Setting.set_IR = 0;
Setting.set_IC = 0;
Setting.values = linspace(Result_ParetoBBSW(2).RelaxSW,...
                        Result_vertex(4).SW - Result_maxrelax(1).SW, round(Setting.Pareto_Point * 0.75));
Setting.primal_set = Result_primal_collect;
tic;
[Result_ParetoSWBB,solution_ParetoSWBB] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoSWBB = toc;
if isfield(Setting,'primal_set')
    Setting = rmfield(Setting, 'primal_set');
end
if isfield(Setting,'Result_primal_Paretostart')
    Setting = rmfield(Setting, 'Result_primal_Paretostart');
end

% j = 0;
% Result_ParetoBBSW_all = [];
% solution_ParetoBBSW_all = [];
% for i = 2:length(Result_ParetoBBSW)
%     j = j+1;
%     Result_ParetoBBSW_all(j) = Result_ParetoBBSW(i);
%     solution_ParetoBBSW_all(j) = solution_ParetoBBSW(i);
% end
j = length(Result_ParetoBBSW);
Result_ParetoBBSW_all = Result_ParetoBBSW(2:j);
solution_ParetoBBSW_all = solution_ParetoBBSW(2:j);
j = j-1;

for i = 2:length(Result_ParetoSWBB)
    j = j+1;
    Result_ParetoBBSW_all(j) = Result_ParetoSWBB(i);
    solution_ParetoBBSW_all(j) = solution_ParetoSWBB(i);
end

plotmatrix = sortrows([[Result_ParetoBBSW_all.welfare]' [Result_ParetoBBSW_all.surplus]'],1);
plot(plotmatrix(:,1), plotmatrix(:,2))
% %% BB_SW_add
% Setting.Pareto_XVar = 'BB';
% Setting.Pareto_YVar = 'SW';
% Setting.set_IR = 0;
% Setting.set_IC = 0;
% Setting.values = linspace(tempsave(1)*1.1-tempsave(2)*0.1,tempsave(2),10);
% 
% tic;
% [Result_ParetoBBSW_add,solution_ParetoBBSW_add] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
% time.ParetoBBSW_add = toc;
% % add
% Result_ParetoBBSW = [Result_ParetoBBSW Result_ParetoBBSW_add(1:11)];
% solution_ParetoBBSW = [solution_ParetoBBSW solution_ParetoBBSW_add(1:11)];
% clear Result_ParetoBBSW_add solution_ParetoBBSW_add

save result_temp
%% BB_IC
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'IC';
Setting.set_IR = 0;
Setting.set_SW = 0;
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);
Setting.Result_primal_Paretostart = Result_maxrelax(1);

tic;
[Result_ParetoBBIC,solution_ParetoBBIC] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBIC = toc;
% 
%% IR_SW
Setting.Pareto_XVar = 'IR';
Setting.Pareto_YVar = 'SW';
Setting.set_BB = 0;
Setting.set_IC = 0;
% Result_maxrelax(2).RelaxIR;
% Result_vertex(5).RelaxIR;   

Setting.values = linspace(Result_maxrelax(2).RelaxIR,Result_vertex(5).RelaxIR,Setting.Pareto_Point);

Setting.Result_primal_Paretostart = Result_maxrelax(2);
tic;
[Result_ParetoIRSW,solution_ParetoIRSW] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoIRSW = toc;


%%
fprintf('基础问题建模总时间%f\n',time.construct);
fprintf('四顶点求解总时间%f\n',time.vertex);
fprintf('最大松弛项求解总时间%f\n',time.maxrelax);
fprintf('BBSW帕累托求解总时间%f\n',time.ParetoBBSW);
fprintf('BBIC帕累托求解总时间%f\n',time.ParetoBBIC);
fprintf('IRSW帕累托求解总时间%f\n',time.ParetoIRSW);
% save_name = [script_type, '_', xlsname ,'_', sheetname, '_', Setting_suffix];
% save(['Result/',save_name], 'Num','Para','Setting','Result_collect','solution_collect','time')
%%
clear Cons_basic Var_basic Cons_num
%%
save_name = [script_type, '_', xlsname ,'_', sheetname, '_', Setting_suffix, '_topology',num2str(Setting.topology)];
save(['Result/',save_name]);
%%
Result_ParetoIRSW = 0;
% save_name = 0;
plot_Pareto_fun(0,solution_ParetoBBSW,Result_ParetoBBSW,Result_ParetoBBIC,Result_ParetoIRSW,save_name)

