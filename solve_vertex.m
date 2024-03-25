function [Result_collect,solution_collect] = solve_vertex(Setting, Para, Num, Var_basic, Cons_basic)
    % 本函数用来求解可行域的四个极点
    Setting.set_IC = 0;
    Setting.set_IR = 0;
    Setting.set_BB = 0;
    Setting.set_SW = 0;
    
    
    
    count = 0;
    fprintf('*****************************验证问题是否可行*********************************')
    % 这里会让所有的都出清, 把各个性质都松弛掉，查看是否可行。
    Setting.choose_relax = 'check_feasible'; %BB,SW,SWmax
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_feasible,solution_feasible] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    count = count+1;
    Result_feasible.name = 'feasible';
    solution_collect(count)= solution_feasible;
    Result_collect(count)= Result_feasible;


%     fprintf('*****************************最大化可能的社会福利*********************************')
%     Setting.choose_relax = 'SWmax'; %BB,SW,SWmax
%     
%     Setting.Iter_init = 1;
%     Setting.Result_init = Setting.Result_primal_SWmax;
% %     Setting.Result_init.RelaxIC = 0;
% %     Setting.Result_init.RelaxIR = 0;
% %     Setting.Result_init.RelaxSW = -Setting.Result_init.SW;
% %     Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
%     
% %     Setting.init_type = 0;
%     
%     [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);      
%     [Result_SWmax,solution_SWmax] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
%     Para.SW_max = Result_SWmax.welfare;%这个离散化了以后会有很多问题
%     count = count+1;
%     Result_SWmax.name = 'SWmax';
%     solution_collect(count)= solution_SWmax;
%     Result_collect(count)= Result_SWmax;
%     % 

    fprintf('*****************************松弛社会福利*********************************')
    Setting.choose_relax = 'SW'; %BB,SW,SWmax
    
    Setting.Iter_init = 1;
    Setting.Result_init = Setting.Result_primal_BBbalance;
%     Setting.Result_init.RelaxIC = 0;
%     Setting.Result_init.RelaxIR = 0;
%     Setting.Result_init.RelaxSW = Result_SWmax.SW - Setting.Result_primal_BBbalance.SW;
%     Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
    
    
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_SW,solution_SW] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    count = count+1;
    Result_SW.name = 'SW';
    solution_collect(count)= solution_SW;
    Result_collect(count)= Result_SW;

    
    %松弛BB,要求达到某个SW.此时的q_ip会往小(绝对值)里计. 我怎么知道
    fprintf('*****************************松弛收支平衡*********************************')
    Setting.choose_relax = 'BB'; %BB,SW,SWmax
    
    Setting.Iter_init = 1;
    Setting.Result_init = Setting.Result_primal_SWmax;
%     Setting.Result_init.RelaxIC = 0;
%     Setting.Result_init.RelaxIR = 0;
%     Setting.Result_init.RelaxSW = 0;
%     Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
    Setting.Result_primal_RelaxBB = Setting.Result_init;
    
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_BB,solution_BB] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    count = count+1;
    Result_BB.name = 'BB';
    solution_collect(count)= solution_BB;
    Result_collect(count)= Result_BB;
    
    
    fprintf('*****************************松弛个体理性*********************************')
    Setting.choose_relax = 'IR'; %BB,SW,SWmax
    if isfield(Setting,'new_version') && Setting.new_version == 1
        Setting.Iter_init = 0;
    else
        Setting.Iter_init = 1;
        Setting.Result_init = Result_BB;%Setting.Result_primal_SWmax;
        Setting.Result_init.RelaxIC = 0;
        Setting.Result_init.RelaxIR = Result_BB.RelaxBB * Para.BBIR_ratio;%Setting.Result_primal_SWmax.RelaxBB * Para.BBIR_ratio;
        Setting.Result_init.RelaxSW = 0;
        Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
        Setting.Result_primal_RelaxIR = Setting.Result_init;
        
    end
    
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_IR,solution_IR] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    count = count+1;
    Result_IR.name = 'IR';
    solution_collect(count)= solution_IR;
    Result_collect(count)= Result_IR;
    
    
    fprintf('*****************************松弛激励相容*********************************')
    Setting.choose_relax = 'IC'; %BB,SW,SWmax
    if isfield(Setting,'new_version') && Setting.new_version == 1
        Setting.Iter_init = 0;
    else
        Setting.Iter_init = 1;
        Setting.Result_init = Result_BB;%Setting.Result_SWmax;
        Setting.Result_init.RelaxIC = Result_BB.RelaxBB * Para.BBIC_ratio;%Setting.Result_SWmax.RelaxBB * Para.BBIC_ratio;
        Setting.Result_init.RelaxIR = 0;
        Setting.Result_init.RelaxSW = 0;
        Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
        Setting.Result_primal_RelaxIC = Setting.Result_init;
    end
    
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_IC,solution_IC] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    count = count+1;
    Result_IC.name = 'IC';
    solution_collect(count)= solution_IC;
    Result_collect(count)= Result_IC;

end