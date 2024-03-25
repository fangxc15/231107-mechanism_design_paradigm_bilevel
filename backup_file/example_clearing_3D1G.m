% 这是一个基于多主体交易的优化问题尝试
clear 

script_type = '3D1G';

Num.I = 4;
Num.P = 20;
Num.C = 4;

Setting_suffix = [num2str(Num.P), 'P'];

Para.SW_max = 0;%这个离散化了以后会有很多问题

Para.type = [1,1,1,2];
Para.qmin = [0,0,0,-1];
Para.qmax = [1,1,1,0];
Para.xmax = [1,1,1,1];
Para.xmin = [0,0,0,0];
Para.distribution = [0,0,0,0];

Para.xpara = [0.6,0.7,0.8,0.2]; %
Para.xstd = [0.1,0.1,0.1,0.1];
Para.company = [1,2,3,4]; %表示分属几号公司.
Para.xworse(Para.type == 1) = Para.xmin(Para.type == 1);
Para.xworse(Para.type == 2) = Para.xmax(Para.type == 2);

Setting.sample_method = 0;
Setting.Pareto_Point = 50;
[Para,Num] = process_Para(Para,Num,Setting);


% % 一般的情况如下,把买方设计在前
% Para.Point = ones(Num.I,1) * [1/(Num.P)*[0:Num.P]];
% % Para.Point = ones(Num.I,1) * [0:Num.I:Num.I*Num.P]/(Num.I*(Num.P+1)-1) +  [0:Num.I-1]'/(Num.I*(Num.P+1)-1) * ones(1,Num.P+1);
% 
% % 如果是起始不同的参数区间
% 
% 
% Para.Interval_len = Para.Point(:,2:Num.P+1) - Para.Point(:,1:Num.P);
% % 这个Probability假设的是均匀分布. 如果从r(x)的积分角度来看，应该采用
% Para.Prob = 0.5 * ([Para.Point(:,2:Num.P+1) 2-Para.Point(:,Num.P+1)] - [-Para.Point(:,1) Para.Point(:,1:Num.P)]); 
% % 买方；卖方
% Para.Dset = find(Para.type == 1); %1是买方
% Para.Gset = find(Para.type == 2); %2是卖方

Para.SW_max = 0;%这个离散化了以后会有很多问题
Setting.Qipp = 1;
Setting.SWtolerance = 1e-5;
Setting.bigM = 10000;

Setting.unified_spread = 1; % 是否设计统一的spread
Setting.riskSW = 1;         % 是否设计带有风险项的SW计算
Setting.topology = 0;
Setting.direct_solve = 0;

%% Basic_construct
tic;
[Var_basic, Cons_basic, Cons_num] = optimization_construct(Num,Para,Setting);
time.construct = toc;
%% Feasible
% [Result_feasible,solution_feasible] = solve_feasible(Setting, Para, Num, Var_basic, Cons_basic, Cons_num);
%% Vertex
tic;
[Result_vertex,solution_vertex] = solve_vertex(Setting, Para, Num, Var_basic, Cons_basic);
time.vertex = toc;
%% maxrelax
tic;
[Result_maxrelax,solution_maxrelax] = solve_maxrelax(Setting, Para, Num, Var_basic, Cons_basic);
time.maxrelax = toc;
%% BB_SW
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'SW';
Setting.set_IR = 0;
Setting.set_IC = 0;
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);
% Setting.values = -1/6:1/120:1/6+1/120; %这是帕累托settign

tic;
[Result_ParetoBBSW,solution_ParetoBBSW] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBSW = toc;
% %% BB_SW_add
% Setting.Pareto_XVar = 'BB';
% Setting.Pareto_YVar = 'SW';
% Setting.set_IR = 0;
% Setting.set_IC = 0;
% Setting.values = [-50/1200:0.1/1200:-49.1/1200  -49/1200:1/1200:-40/1200];
% 
% tic;
% [Result_ParetoBBSW_add,solution_ParetoBBSW_add] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
% time.ParetoBBSW_add = toc;
% %% add
% Result_ParetoBBSW = [Result_ParetoBBSW Result_ParetoBBSW_add(1:12)];
% solution_ParetoBBSW = [solution_ParetoBBSW solution_ParetoBBSW_add(1:12)];
% clear Result_ParetoBBSW_add solution_ParetoBBSW_add
%% BB_IC
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'IC';
Setting.set_IR = 0;
Setting.set_SW = 0;
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);

% Setting.values = 0-1/60:1/60:1/6+1/60;

tic;
[Result_ParetoBBIC,solution_ParetoBBIC] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBIC = toc;

%% IR_SW
Setting.Pareto_XVar = 'IR';
Setting.Pareto_YVar = 'SW';
Setting.set_BB = 0;
Setting.set_IC = 0;
Setting.values = linspace(Result_maxrelax(2).RelaxIR,Result_vertex(5).RelaxIR,Setting.Pareto_Point);

tic;
[Result_ParetoIRSW,solution_ParetoIRSW] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoIRSW = toc;

%% 输出时间
fprintf('基础问题建模总时间%f\n',time.construct);
fprintf('四顶点求解总时间%f\n',time.vertex);
fprintf('最大松弛项求解总时间%f\n',time.maxrelax);
fprintf('BBSW帕累托求解总时间%f\n',time.ParetoBBSW);
fprintf('BBIC帕累托求解总时间%f\n',time.ParetoBBIC);
fprintf('IRSW帕累托求解总时间%f\n',time.ParetoIRSW);

%% 保存
clear Cons_basic Var_basic Cons_num
save_name = [script_type,'_', Setting_suffix];
save(['Result/',save_name]);
%% 画图
plot_Pareto_fun(1,solution_ParetoBBSW,Result_ParetoBBSW,Result_ParetoBBIC,Result_ParetoIRSW,save_name)


