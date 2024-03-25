clear 
script_type = 'xls';
Setting.new_version = 2;
version_suffix = '_newversion2';
folder_name = ['Result',version_suffix];
mkdir(folder_name);


%%
xlsname = 'Clearing_IEEE30';
sheetname = 'Version2'; %这个地方开始出现不平衡的原因是把抽样方式改了。采用了间隔掺入式的抽样
IEEE_30_Setting;

%% Heursitic find feasibile (find max SW, s.t. BB= 0) Optimization primal 
Result_maxSW = optimization_primal(Num,Para,Setting, zeros(1,Num.I));
Para.SW_max = Result_maxSW.SW;
Result_primal_BBbalance = heuristic_find_BBbalance(Setting.heuristic_spread_start,...
    Setting.heuristic_spread_end, Setting.heuristic_spread_tolerance,Num,Para,Setting);

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
Result_primal_BBmax = heuristic_find_BBmax(Setting.heuristic_spread_start,...
    Setting.heuristic_spread_end, Setting.heuristic_spread_tolerance,Num,Para,Setting);

BBmaxno = find(([Result_primal_BBmax.BB] == max([Result_primal_BBmax.BB])));
Setting.Result_primal_BBmax = Result_primal_BBmax(BBmaxno(1));


%% Basic
tic;
[Var_basic, Cons_basic] = optimization_construct(Num,Para,Setting);
time.construct = toc;
% %% Optimization primal
% Setting.Result_primal = optimization_primal(Num,Para,Setting);

%% Vertex
tic;
Setting.hard_assign = [0]; %这个表示不要直接assign给他很粗暴的值

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
% Setting.Result_primal.SW


%% maxrelax
tic;
[Result_maxrelax,solution_maxrelax] = solve_maxrelax(Setting, Para, Num, Var_basic, Cons_basic);
time.maxrelax = toc;
%% More heursitic solutions
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

%% 1 BB_SW
% 在IR,IC = 0的情况下，之前求解过，哪些部分是BB-SW的Pareto前沿？这个值得考虑
% 分别是RelaxBB最小的点和RelaxSW最小的点。这两个点有意义。
%
Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'SW';
Setting.set_IR = 0;
Setting.set_IC = 0;
% Result_maxrelax(1).RelaxBB; %代表最大的可能盈余
% Result_vertex(3).RelaxBB;   %代表要求保证其他项时的盈余relax项
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(3).RelaxBB,Setting.Pareto_Point);
tic;
[Result_ParetoBBSW,solution_ParetoBBSW] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBSW = toc;

plotmatrix = sortrows([[Result_ParetoBBSW.welfare]' [Result_ParetoBBSW.surplus]' [Result_ParetoBBSW.kIR]'],1);
% plotmatrix = plotmatrix(find([solution_ParetoBBSW.problem] == 0),:);
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
plot(plotmatrix(:,1), plotmatrix(:,2))
%% 2 BB_IR
% 在SW,IC = 0的情况下，之前求解过，哪些部分是BB-IR的Pareto前沿？这个值得考虑
% 分别是Relax_BB最小的点和Relax_IR最小的点。这两个点有意义。
% 但其实也可以存在这样一个点，让Relax_BB更小，Relax_IR更大。

Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'IR';
Setting.set_IC = 0;
Setting.set_SW = 0;
% Result_vertex(4).RelaxBB Result_vertex(4).RelaxIR  %0    0.1856
% Result_vertex(3).RelaxBB Result_vertex(3).RelaxIR  %5450 0
% 其实最好的选择点是这样的

%maxrelax(1).RelaxBB根本不是在SW=0的时候达成的
Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(3).RelaxBB,Setting.Pareto_Point);

tic;
[Result_ParetoBBIR,solution_ParetoBBIR] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBIR = toc;

plotmatrix = sortrows([[Result_ParetoBBIR.RelaxBB]' [Result_ParetoBBIR.RelaxIR]' ...
    [Result_ParetoBBIR.kIR]'],1);
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
%为什么会出现直线问题? 因为本身的出清结果是已经固定的。由于边际上的RX没有冲抵的区间，所以只能这样。
%修改IR的定义方式后解决了该问题
plot(plotmatrix(:,1), plotmatrix(:,2))
% 感觉是简直可以无限下去
%% 3 BB_IC
% 在SW,IR = 0的情况下，之前求解过，哪些部分是BB-IC的Pareto前沿？这个值得考虑
% 分别是Relax_BB最小的点和Relax_IC最小的点。这两个点有意义。
% 但其实也可以存在这样一个点，让Relax_BB更小，Relax_IR更大。

Setting.Pareto_XVar = 'BB';
Setting.Pareto_YVar = 'IC';
Setting.set_IR = 0;
Setting.set_SW = 0;
% Result_vertex(5).RelaxBB Result_vertex(5).RelaxIC  %infeasible
% Result_vertex(3).RelaxBB Result_vertex(3).RelaxIC  %5450 0
% 其实最好的选择点是这样的

Setting.values = linspace(Result_maxrelax(1).RelaxBB,Result_vertex(3).RelaxBB,Setting.Pareto_Point);

tic;
[Result_ParetoBBIC,solution_ParetoBBIC] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoBBIC = toc;

plotmatrix = [[Result_ParetoBBIC.RelaxBB]' [Result_ParetoBBIC.RelaxIC]' [Result_ParetoBBIC.kIC]'];
% plotmatrix = plotmatrix(find([solution_ParetoBBIC.problem] == 0),:);
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
plotmatrix = sortrows(plotmatrix,1);
%在RelaxBB>0的时候，此时也能够Relax部分IC, 因为此时也能让k^IR有部分补偿
%这里应该要把部分infeasible去除掉
plot(plotmatrix(:,1), plotmatrix(:,2))

%% 4 SW-IR
Setting.Pareto_XVar = 'SW';
Setting.Pareto_YVar = 'IR';
Setting.set_BB = 0;
Setting.set_IC = 0;

% 在BB,IC = 0的情况下，之前求解过，哪些部分是SW-IR的Pareto前沿？这个值得考虑
% 分别是Relax_SW最小的点和Relax_IR最小的点。这两个点有意义。
% Result_vertex(2).RelaxIR Result_vertex(2).RelaxSW  %0         158.1263
% Result_vertex(4).RelaxIR Result_vertex(4).RelaxSW  %0.1856    0
% 其实最好的选择点是这样的

% 这是SW-IR的Setting.Point
% Setting.values = linspace(Result_maxrelax(2).RelaxIR,Result_vertex(4).RelaxIR,Setting.Pareto_Point);
Setting.values = linspace(0,Result_vertex(2).RelaxSW * 2,Setting.Pareto_Point); 
%实际上根本不是这样，在RelaxSW后肯定IR始终都是0

tic;
[Result_ParetoSWIR,solution_ParetoSWIR] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoSWIR = toc;
% 

plotmatrix = [[Result_ParetoSWIR.RelaxSW]' [Result_ParetoSWIR.RelaxIR]' [Result_ParetoSWIR.kIR]'];
% plotmatrix = plotmatrix(find([solution_ParetoSWIR.problem] == 0),:);
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
plotmatrix = sortrows(plotmatrix,1);
%在RelaxBB>0的时候，此时也能够Relax部分IC, 因为此时也能让k^IR有部分补偿
%这里应该要把部分infeasible去除掉
plot(plotmatrix(:,1), plotmatrix(:,2))



%% 5 SW-IC
Setting.Pareto_XVar = 'SW';
Setting.Pareto_YVar = 'IC';
Setting.set_BB = 0;
Setting.set_IR = 0;

% 在BB,SW = 0的情况下，之前求解过，哪些部分是IR-IC的Pareto前沿？这个值得考虑
% 是RelaxIR最大的点和RelaxIR为0的点
% Result_vertex(4).RelaxIR
% 其实最好的选择点是这样的


Setting.values = linspace(0,Result_vertex(2).RelaxSW * 2,Setting.Pareto_Point);


tic;
[Result_ParetoSWIC,solution_ParetoSWIC] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoSWIC = toc;
% 

plotmatrix = [[Result_ParetoSWIC.RelaxSW]' [Result_ParetoSWIC.RelaxIC]'];
% plotmatrix = plotmatrix(find([solution_ParetoSWIC.problem] == 0),:);
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);

plotmatrix = sortrows(plotmatrix,1);
%在RelaxBB>0的时候，此时也能够Relax部分IC, 因为此时也能让k^IR有部分补偿
%这里应该要把部分infeasible去除掉
plot(plotmatrix(:,1), plotmatrix(:,2))
%% 6 IR-IC
Setting.Pareto_XVar = 'IR';
Setting.Pareto_YVar = 'IC';
Setting.set_BB = 0;
Setting.set_SW = 0;

% 在BB,SW = 0的情况下，之前求解过，哪些部分是IR-IC的Pareto前沿？这个值得考虑
% 是RelaxIR最大的点和RelaxIR为0的点
% Result_vertex(4).RelaxIR
% 其实最好的选择点是这样的


Setting.values = linspace(Result_maxrelax(2).RelaxIR,Result_vertex(4).RelaxIR,Setting.Pareto_Point);


tic;
[Result_ParetoIRIC,solution_ParetoIRIC] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
time.ParetoIRIC = toc;
% 

plotmatrix = [[Result_ParetoIRIC.RelaxIR]' [Result_ParetoIRIC.RelaxIC]' ...
                [Result_ParetoIRIC.kIR]' [Result_ParetoIRIC.kIC]'];
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0,:);
% plotmatrix = plotmatrix(find([solution_ParetoIRIC.problem] == 0),:);
plotmatrix = sortrows(plotmatrix,1);
%在RelaxBB>0的时候，此时也能够Relax部分IC, 因为此时也能让k^IR有部分补偿
%这里应该要把部分infeasible去除掉
plot(plotmatrix(:,1), plotmatrix(:,2))



%% 求SW,BB,IR三者的Pareto Frontier
matrix_BBSW = [[Result_ParetoBBSW.RelaxSW]' [Result_ParetoBBSW.RelaxBB]' (1:length(Result_ParetoBBSW))'];
matrix_BBSW = matrix_BBSW(matrix_BBSW(:,1) ~= 0 | matrix_BBSW(:,2) ~= 0,:);
matrix_BBSW = matrix_BBSW(matrix_BBSW(:,2)>=0,:);
% matrix_BBSW = [Result_vertex(2).RelaxSW 0;matrix_BBSW];

%固定SW的时候进行BB/IR的权衡
time.ParetoSWIRBB = 0;
Result_ParetoSWIRBB = [];
solution_ParetoSWIRBB = [];
if isfield(Setting,'primal_set')
    Setting = rmfield(Setting,'primal_set');
end
for i = 1:length(matrix_BBSW)
    Setting.Pareto_XVar = 'IR';
    Setting.Pareto_YVar = 'BB';
    Setting.set_IC = 0;
    Setting.set_SW = matrix_BBSW(i,1);
    
    Setting.Result_primal_Paretostart = Result_ParetoBBSW(matrix_BBSW(i,3));
    
    Setting.values = linspace(0,Result_vertex(4).RelaxIR,10);
    
    tic;
    [Result_temp,solution_temp] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic);
    time.ParetoSWIRBB = time.ParetoSWIRBB + toc;
    k = length(Result_ParetoSWIRBB);
    if k == 0
        Result_ParetoSWIRBB = Result_temp;
        solution_ParetoSWIRBB = solution_temp;
    else
        for j = 1:length(Result_temp)
            k = k+1;
            Result_ParetoSWIRBB(k) = Result_temp(j);
            solution_ParetoSWIRBB(k) = solution_temp(j);
        end
    end  
end

plotmatrix = [[Result_ParetoSWIRBB.RelaxSW]' [Result_ParetoSWIRBB.RelaxIR]' ...
                [Result_ParetoSWIRBB.RelaxBB]'];
% plotmatrix = plotmatrix(find([solution_ParetoSWIRBB.problem] == 0),:);
plotmatrix = plotmatrix(plotmatrix(:,1) ~= 0 | plotmatrix(:,2) ~= 0 | ...
    plotmatrix(:,3) ~= 0 ,:);

plotmatrix(:,1) = plotmatrix(:,1)/Para.SW_max;
plotmatrix(:,3) = plotmatrix(:,3)/Para.SW_max;

set(groot,'defaultfigurePosition',[300 300 460 390]);
set(0,'defaultfigurecolor','w'); %设置背景颜色为白色
set(groot,'defaultAxesFontWeight','bold');
set(groot,'defaultLegendFontSize',12);
set(groot,'defaultAxesFontSize',12);

[X,Y] = meshgrid(linspace(min(plotmatrix(:,1)),max(plotmatrix(:,1)),20), ...
        linspace(min(plotmatrix(:,2)),max(plotmatrix(:,2)),20));
Z = griddata(plotmatrix(:,1),plotmatrix(:,2),plotmatrix(:,3),X,Y);
s = surf(X,Y,Z);
s.FaceAlpha = 0.9;
s.EdgeColor = 'none';
s.FaceColor = 'interp';
set(gca,'XDir','reverse')
set(gca,'YDir','reverse')
xlabel('社会福利松弛');
% ylim([0,0.08]);
ylabel('个体理性松弛');
zlabel('收支平衡松弛');
title('帕累托前沿')

Picture_root_folder = ['Picture', version_suffix];  
mkdir(Picture_root_folder);
Picture_folder = [Picture_root_folder,'/',save_name];
mkdir(Picture_folder);

saveas(1,[Picture_folder,'/','SWIRBB_Pareto.jpg'])

% meshX = reshape(plotmatrix(:,1),10,length(matrix_BBSW));
% meshY = reshape(plotmatrix(:,2),10,length(matrix_BBSW));
% meshZ = reshape(plotmatrix(:,3),10,length(matrix_BBSW));
% % meshZ(meshZ<0) = 0;
% s = surf(meshX,meshY,meshZ);
% % axis equal
% % axis tight
% s.FaceAlpha = 0.9;
% s.EdgeColor = 'none';
% s.FaceColor = 'interp';
%%
fprintf('基础问题建模总时间%f\n',time.construct);
fprintf('四顶点求解总时间%f\n',time.vertex);
fprintf('最大松弛项求解总时间%f\n',time.maxrelax);
fprintf('BBSW帕累托求解总时间%f\n',time.ParetoBBSW);
fprintf('BBIC帕累托求解总时间%f\n',time.ParetoBBIC);
fprintf('BBIR帕累托求解总时间%f\n',time.ParetoBBIR);
fprintf('IRSW帕累托求解总时间%f\n',time.ParetoIRSW);
fprintf('IRIC帕累托求解总时间%f\n',time.ParetoIRIC);
fprintf('SWIC帕累托求解总时间%f\n',time.ParetoSWIC);

% save_name = [script_type, '_', xlsname ,'_', sheetname, '_', Setting_suffix];
% save(['Result/',save_name], 'Num','Para','Setting','Result_collect','solution_collect','time')
%%
clear Cons_basic Var_basic Cons_num
save_name = [script_type, '_', xlsname ,'_', sheetname, '_', Setting_suffix, '_topology',num2str(Setting.topology)];
save([folder_name ,'/',save_name]);
%%
% plot_Pareto_fun(1,solution_ParetoBBSW,Result_ParetoBBSW,Result_ParetoBBIC,Result_ParetoIRSW,save_name)

