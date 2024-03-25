clear 
script_type = 'xls';
xlsname = 'Clearing_IEEE30';
sheetname = 'Version1'; %这个地方开始出现不平衡的原因是把抽样方式改了。采用了间隔掺入式的抽样

Setting.NumP = 5;
Setting_suffix = [num2str(Setting.NumP), 'P'];
Setting.sample_method = 2; %0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
Setting.cumulative_curve = 1;
Setting.percent_strategic = 0.2;


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
    Para.BranchMax = Para.BranchMax(sum(abs(Para.GSDF),2) > 0);
    Para.GSDF = Para.GSDF(sum(abs(Para.GSDF),2) > 0,:);
%     Para.BranchMax(find(Para.BranchMax == 0)) = 1000;
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

% 
% % temp_D = Para.Point - min(Para.Point,[],2);
% % temp_G = max(Para.Point,[],2) -Para.Point;
% temp_D = Para.Point - Para.xmin';
% temp_G = Para.xmax' - Para.Point;
% 
% 
% 
% % 
% % % 是平衡的
% % sum(sum(Result_vertex(5).mip .* Para.Prob))
% % sum(sum(Result_vertex(6).mip .* Para.Prob))
% % 
% % % 概率总额
% % sum(Para.Prob,2)
% % 
% % Result_vertex(4).RelaxBB
% % Result_vertex(5).RelaxIR
% % Result_vertex(6).RelaxIC
% % 
% % Result_vertex(4).rip %BB
% % Result_vertex(5).rip %IR
% % Result_vertex(6).rip %IC
% % 
% % sum(Result_vertex(4).rip .* Para.Prob,2)
% % sum(Result_vertex(5).rip .* Para.Prob,2)
% % sum(Result_vertex(4).rip .* Para.Prob,2) - sum(Result_vertex(5).rip .* Para.Prob,2)
% % 
% % sum(sum(Result_vertex(4).rip .* Para.Prob))
% % sum(sum(Result_vertex(5).rip .* Para.Prob))
% % sum(sum(Result_vertex(6).rip .* Para.Prob))
% % % 
% % 
% % (sum(sum(Result_vertex(4).rip .* Para.Prob)) - sum(sum(Result_vertex(5).rip .* Para.Prob)))/Num.I
% % 
% % % 平均值就是0.1297;
% % (Result_vertex(4).rip - Result_vertex(6).rip)./[temp_D(Para.Dset,:);temp_G(Para.Gset,:)]
% % totalIC_temp = sum(sum([temp_D(Para.Dset,:);temp_G(Para.Gset,:)] .* Para.Prob,2)) % 就是因为这个上下不一样
% % (sum(sum(Result_vertex(4).rip .* Para.Prob)) - sum(sum(Result_vertex(6).rip .* Para.Prob)))/totalIC_temp
% % 
% %% enumerate_spread
% 
% [result_collect_spread,new_Para] = enumerate_spread(Para,Num,[0:0.2:10],1); %分别是spread_vector和resolution(精度)
% % for i = 1:Num.I
% %     sum(((new_Para.Point(i,:) - new_Para.xmin(i)) .* new_Para.Prob(i,:)))
% % end
% new_temp_D = new_Para.Point - new_Para.xmin';
% new_temp_G = new_Para.xmax' - new_Para.Point;
% % new_Para抽样是均匀的，所以可以这么做
% new_Para.IRIC_ratio = Num.I/sum(sum([new_temp_D(Para.Dset,:);new_temp_G(Para.Gset,:)] .* new_Para.Prob,2)); %这个比较接近实值.是一个比较好的状态;
% new_Para.BBIR_ratio = 1/Num.I;
% new_Para.BBIC_ratio = 1/sum(sum([new_temp_D(Para.Dset,:);new_temp_G(Para.Gset,:)] .* new_Para.Prob,2));


%% Heursitic find feasibile (find max SW, s.t. BB= 0) Optimization primal 
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
tic;
% Setting.Result_BBbalance = Result_vertex(3);
[Result_maxrelax,solution_maxrelax] = solve_maxrelax(Setting, Para, Num, Var_basic, Cons_basic);
time.maxrelax = toc;

for i = 1:length(Result_maxrelax)
    Result_maxrelax(i).RelaxSW = Result_vertex(4).SW - Result_maxrelax(i).SW;
end

%% BB_SW
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'SW';
Setting.set_IR = 0;
Setting.set_IC = 0;
% Result_maxrelax(1).RelaxBB; %代表最大的可能盈余
% Result_vertex(4).RelaxBB;   %代表要求保证其他项时的盈余relax项
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);
% tempsave = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,30);
Setting.Result_primal_Paretostart = Result_maxrelax(1);

tic;
[Result_ParetoBBSW,solution_ParetoBBSW] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBSW = toc;
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
% 

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
save_name = [script_type, '_', xlsname ,'_', sheetname, '_', Setting_suffix, '_topology',num2str(Setting.topology)];
save(['Result/',save_name]);
%%
plot_Pareto_fun(1,solution_ParetoBBSW,Result_ParetoBBSW,Result_ParetoBBIC,Result_ParetoIRSW,save_name)

