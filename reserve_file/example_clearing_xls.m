clear 
script_type = 'xls';
xlsname = 'Clearing_V1_node';
sheetname = '8asymmetry_monopoly'; %这个地方开始出现不平衡的原因是把抽样方式改了。采用了间隔掺入式的抽样

Setting.NumP = 10;
Setting_suffix = [num2str(Setting.NumP), 'P'];
Setting.sample_method = 1; %0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
Setting.cumulative_curve = 1;
[Para,Num] = read_info([xlsname,'.xlsx'],sheetname, Setting);
%% 考虑拓扑

Setting.topology = 0;
if Setting.topology == 1
    %导入拓扑的程序
    mpc = case5;
    GSDF_lb = makePTDF(mpc); %node和Branch之间的转换矩阵, Num.Branch *Num.Node.
    Para.GSDF = GSDF_lb(:,Para.node); %主体和Branch之间的转换矩阵, Num.Branch *Num.Participant.
    Para.BranchMax = mpc.branch(:,6);
    Para.BranchMax = ones(size(Para.BranchMax,1),1);
    Num.Branch = size(mpc.branch,1);
end
%%
Para.SW_max = 0;%这个离散化了以后会有很多问题
Setting.Qipp = 0;
Setting.SWtolerance = 1e-5;
Setting.Complementary_bigM = max(abs(min(Para.qmin)),abs(max(Para.qmax)));
Setting.Complementary_bigM = max(max(Para.xmax) * 2, Setting.Complementary_bigM);
if Setting.topology == 1
    Setting.Complementary_bigM = max(Setting.Complementary_bigM, max(Para.BranchMax) * 2)*2;
end

Setting.bigM = 10000;

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


%% Basic
tic;
[Var_basic, Cons_basic] = optimization_construct(Num,Para,Setting);
time.construct = toc;
%% Vertex
tic;
[Result_vertex,solution_vertex] = solve_vertex(Setting, Para, Num, Var_basic, Cons_basic);
time.vertex = toc;

Para.BBIR_ratio = 1/Num.C; %Result_vertex(5).RelaxIR/Result_vertex(4).RelaxBB

for c = 1:Num.C
    tempno = Para.company_set(c).no;
    temp_len(c) = sum(sum(Para.Interval_len(tempno,:) .^2,1).^0.5);
end         
Para.BBIC_ratio = 1/sum(temp_len); %Result_vertex(6).RelaxIC/Result_vertex(4).RelaxBB

%% maxrelax
tic;
[Result_maxrelax,solution_maxrelax] = solve_maxrelax(Setting, Para, Num, Var_basic, Cons_basic);
time.maxrelax = toc;
%% BB_SW
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'SW';
Setting.set_IR = 0;
Setting.set_IC = 0;
% Result_maxrelax(1).RelaxBB; %代表最大的可能盈余
% Result_vertex(4).RelaxBB;   %代表要求保证其他项时的盈余relax项
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,Setting.Pareto_Point);
% tempsave = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(4).RelaxBB,30);
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

tic;
[Result_ParetoBBIC,solution_ParetoBBIC] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBIC = toc;
% 
% IR_SW
Setting.Pareto_XVar = 'IR';
Setting.Pareto_YVar = 'SW';
Setting.set_BB = 0;
Setting.set_IC = 0;
% Result_maxrelax(2).RelaxIR;
% Result_vertex(5).RelaxIR;   

Setting.values = linspace(Result_maxrelax(2).RelaxIR,Result_vertex(5).RelaxIR,Setting.Pareto_Point);


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

