function [Result_collect,solution_collect] = solve_feasible(Setting, Para, Num, Var_basic, Cons_basic, Cons_num)

    count = 0;
    fprintf('*****************************验证问题是否可行*********************************')
    % 这里会让所有的都出清, 把各个性质都松弛掉，查看是否可行。
%     Setting.choose_relax = 'check_feasible'; %BB,SW,SWmax
%     [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    
    Var = Var_basic;
    Cons = Cons_basic;
    
    % 这里是我们提出的check_feasible的内容.按照我们所提的方法进行出清
    Cons(Cons_num.qbalance) = [];
    Cons(reshape(Cons_num.Qqlink,1,numel(Cons_num.Qqlink))) = [];
    
    Cons = [Cons, Var.qip(1,:) == 32/27*max(0,Para.Point(1,:) - 0.25)]; %这个地方写的不对
    Cons = [Cons, Var.qip(2,:) == -max(0,0.75 - 27/32 * Para.Point(2,:))];
    obj = 0;

    [Result_feasible,solution_feasible] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    count = count+1;
    Result_feasible.name = 'feasible';
    solution_collect(count)= solution_feasible;
    Result_collect(count)= Result_feasible;