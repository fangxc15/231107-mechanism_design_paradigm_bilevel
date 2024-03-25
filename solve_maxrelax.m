function [Result_collect,solution_collect] = solve_maxrelax(Setting, Para, Num, Var_basic, Cons_basic)
    
    Setting.set_IC = 0;
    Setting.set_IR = 0;
    Setting.set_BB = 0;
    Setting.set_SW = 0;
    count = 0;
    
    
    fprintf('*****************************最大化可能的盈余项*********************************')
    Setting.choose_relax = 'BBmax'; %BB,SW,SWmax
    
    Setting.Iter_init = 1;
    Setting.Result_init = Setting.Result_primal_BBmax; %Setting.Result_primal;
%     Setting.Result_init.RelaxIC = 0;
%     Setting.Result_init.RelaxIR = 0;
%     Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
%     Setting.BBmax = Setting.Result_init;
    
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_BBmax,solution_BBmax] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
   
    count = count+1;
    Result_BBmax.name = 'BBmax';
    solution_collect(count)= solution_BBmax;
    Result_collect(count)= Result_BBmax;
    
    fprintf('*****************************最大化IR松弛项*********************************')
    Setting.choose_relax = 'IRmax'; %BB,SW,SWmax
    
    if isfield(Setting,'new_version') && Setting.new_version == 1
        Setting.Iter_init = 0;
    else  
        Setting.Iter_init = 1;
        Setting.Result_init = Result_BBmax;%Setting.Result_primal_BBmax; 
        Setting.Result_init.RelaxIC = 0;
        Setting.Result_init.RelaxIR = Result_BBmax.RelaxBB * Para.BBIR_ratio;%Setting.Result_BBmax.RelaxBB * Para.BBIR_ratio;
        Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
        Setting.Result_primal_IRmax = Setting.Result_init;
    end
        
    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_IRmax,solution_IRmax] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    
    count = count+1;
    Result_IRmax.name = 'IRmax';
    solution_collect(count)= solution_IRmax;
    Result_collect(count)= Result_IRmax;
    
    fprintf('*****************************最大化IC松弛项*********************************')
    Setting.choose_relax = 'ICmax'; %BB,SW,SWmax
    if isfield(Setting,'new_version') && Setting.new_version == 1
        Setting.Iter_init = 0;
    else
        Setting.Iter_init = 1;
        Setting.Result_init = Result_BBmax;%Setting.Result_primal_BBmax; 
        Setting.Result_init.RelaxIC = Result_BBmax.RelaxBB * Para.BBIC_ratio; %Setting.Result_BBmax.RelaxBB * Para.BBIC_ratio;
        Setting.Result_init.RelaxIR = 0;
        Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
        Setting.Result_primal_ICmax = Setting.Result_init;
    end

    [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
    [Result_ICmax,solution_ICmax] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
    
    count = count+1;
    Result_ICmax.name = 'ICmax';
    solution_collect(count)= solution_ICmax;
    Result_collect(count)= Result_ICmax;
    
end