function [Result_collect,solution_collect] = solve_Pareto(Setting, Para, Num, Var_basic, Cons_basic)
    % 求解帕累托前沿的方法, 需要明确求解哪个帕累托
    
    count = 0;
%     fprintf('*****************************最大化可能的社会福利*********************************')
%     Setting.choose_relax = 'SWmax'; %BB,SW,SWmax
%     Setting.Iter_init = 1;
%     Setting.Result_init = Setting.Result_primal_SWmax;
% %     Setting.Result_init.RelaxIC = 0;
% %     Setting.Result_init.RelaxIR = 0;
% %     Setting.Result_init.RelaxSW = -Setting.Result_init.SW;
% %     Setting.Result_init = calculate_MX(Setting.Result_init, Setting, Num, Para);
%     
%     
%     [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
%     [Result_SWmax,solution_SWmax] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
%     Para.SW_max = Result_SWmax.welfare;%这个离散化了以后会有很多问题 
%     Result_SWmax.name = 'SWmax';
%     index = index + 1;
%     solution_collect(count)= solution_SWmax;
%     Result_collect(count)= Result_SWmax;
%     Setting = rmfield(Setting,'Result_init');
    
    
    Setting.choose_relax = Setting.Pareto_YVar;
    
    if ~isfield(Setting,'set_SW')
        Setting.set_SW = 0;
    end
    if ~isfield(Setting,'set_BB')
        Setting.set_BB = 0;
    end
    if ~isfield(Setting,'set_IR')
        Setting.set_IR = 0;
    end
    if ~isfield(Setting,'set_IC')
        Setting.set_IC = 0;
    end
    
    if isfield(Setting, 'primal_set')
        for j = 1:length(Setting.primal_set)
            Setting.primal_set(j).RelaxSW = Para.SW_max - Setting.primal_set(j).SW;
        end

        if strcmp(Setting.Pareto_XVar , 'BB')
            candidate = [Setting.primal_set.RelaxBB];
        elseif strcmp(Setting.Pareto_XVar , 'SW')
            candidate = [Setting.primal_set.RelaxSW];
        elseif strcmp(Setting.Pareto_XVar , 'IC')
            candidate = [Setting.primal_set.RelaxIC];
        elseif strcmp(Setting.Pareto_XVar , 'IR')
            candidate = [Setting.primal_set.RelaxIR];
        end

        candidate = [candidate' (1:length(candidate))'];
        candidate = sortrows(candidate,1,'descend');
    end

     
    for i = 1:length(Setting.values)
        fprintf('*****************************%s=%f,求%s的帕累托前沿*********************************',Setting.Pareto_XVar,Setting.values(i),Setting.Pareto_YVar)
        Setting = setfield(Setting,['set_',Setting.Pareto_XVar],Setting.values(i));
        [Var, Cons, obj] = optimization_tune(Setting,Para, Var_basic, Cons_basic);
        
        Setting.Iter_init = 1;
        if i == 1    
            if isfield(Setting,'Result_primal_Paretostart')
                Setting.Result_init = Setting.Result_primal_Paretostart;
                Setting.Result_init.RelaxSW = Para.SW_max - Setting.Result_init.SW;
            end
        else
            Setting.Result_init = Result_temp;
        end
        
        if isfield(Setting, 'primal_set')
            choose_num = min(find(candidate(:,1) <= Setting.values(i)));
            if ~isempty(choose_num)
                if isfield(Setting,'Result_init')
                    if candidate(choose_num,1) > getfield(Setting.Result_init,['Relax',Setting.Pareto_XVar]) %如果比他更紧
                        Setting.Result_init = Setting.primal_set(candidate(choose_num,2));
                    end
                else
                    Setting.Result_init = Setting.primal_set(candidate(choose_num,2));
                end
            end
        end
        
        [Result_temp,solution_temp] = optimization_solve(Num,Para,Setting,Var, Cons, obj);
        Result_temp.name = [Setting.Pareto_XVar,' = ',num2str(Setting.values(i),4),' min ',Setting.Pareto_YVar];
        count = count + 1;
        Result_collect(count)= Result_temp;
        solution_collect(count)= solution_temp;
    end 
        
 
end