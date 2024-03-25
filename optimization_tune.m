function  [Var, Cons, obj] = optimization_tune(Setting, Para, Var_basic, Cons_basic)

    Var = Var_basic;
    Cons = Cons_basic;
%     if strcmp(Setting.choose_relax,'Pareto')
%         % 如果采用帕累托的方式,几个帕累托
%         
%         if Setting.
%         if isfield(Setting,'set_SW')
%             if ~isnan(Setting.set_SW)
%                 Cons = [Cons, Var.RelaxSW == Setting.set_SW];
%             end 
%         end 
%         
%         if isfield(Setting,'set_BB')
%             if ~isnan(Setting.set_BB)
%                 Cons = [Cons, Var.RelaxBB == Setting.set_BB];
%             end 
%         end 
        
        
            
    if strcmp(Setting.choose_relax,'check_feasible')
%         Cons = [Cons, Var.Qipp(:,601) == [1;1;-1;-1]];
        % 这里的feasible相当于全部出清; 但是如果买方卖方的个数不等，全部出清本来也无法保证收支平衡吧
%         if isfield(Var,'Qipp')
%             Cons = [Cons, Var.Qipp(Para.Gset, :) == Para.qmin(Para.Gset)' * ones(1,size(Var.Qipp,2))];
%             Cons = [Cons, Var.Qipp(Para.Dset, :) == Para.qmax(Para.Dset)' * ones(1,size(Var.Qipp,2))];
% %         else
%             Cons = [Cons, Var.qip(Para.Gset, :) == Para.qmin(Para.Gset)' * ones(1,size(Var.qip,2))];
%             Cons = [Cons, Var.qip(Para.Dset, :) == Para.qmax(Para.Dset)' * ones(1,size(Var.qip,2))];                   
%         end
        obj = 0;% 应该让它全部不出请才对
        Cons = [Cons, Var.QX(Para.Gset, :) == Para.qmax(Para.Gset)' * ones(1,size(Var.QX,2))];
        Cons = [Cons, Var.QX(Para.Dset, :) == Para.qmin(Para.Dset)' * ones(1,size(Var.QX,2))];       

%         obj = Var.RelaxSW + 10000 * (Var.RelaxIC + Var.RelaxIR);

    elseif strcmp(Setting.choose_relax,'SW')
        Cons = [Cons, Var.RelaxBB <= Setting.set_BB, Var.RelaxIC <= Setting.set_IC, Var.RelaxIR <= Setting.set_IR];
        obj = Var.RelaxSW;
        
    elseif strcmp(Setting.choose_relax,'BB')
        Cons = [Cons, Var.RelaxSW <= Setting.set_SW, Var.RelaxIC <= Setting.set_IC, Var.RelaxIR <= Setting.set_IR];
        obj = Var.RelaxBB;
        
    elseif strcmp(Setting.choose_relax,'IC')
        Cons = [Cons, Var.RelaxSW <= Setting.set_SW, Var.RelaxBB <= Setting.set_BB, Var.RelaxIR <= Setting.set_IR];
        obj = Var.RelaxIC;
        
    elseif strcmp(Setting.choose_relax,'IR')
        Cons = [Cons, Var.RelaxSW <= Setting.set_SW, Var.RelaxBB <= Setting.set_BB, Var.RelaxIC <= Setting.set_IC];
        obj = Var.RelaxIR;
    
        
    elseif strcmp(Setting.choose_relax,'SWmax') % 指的是允许收支不平衡，尝试最大化社会福利,本来的边际出清思路
        Cons = [Cons, Var.RelaxIC <= Setting.set_IC, Var.RelaxIR <= Setting.set_IR];
        obj = Var.RelaxSW; %最小化松弛量，相当于最大化社会福利
        Cons = [Cons, Var.RelaxSW >= -Setting.bigM];
        
    elseif strcmp(Setting.choose_relax,'BBmax') % 指的是允许社会福利不够，尝试最大化BB的盈余
        Cons = [Cons, Var.RelaxIC <= Setting.set_IC, Var.RelaxIR <= Setting.set_IR];
        obj = Var.RelaxBB; %最小化松弛量，相当于最大化盈余量
        Cons = [Cons, Var.RelaxBB >= -Setting.bigM];
        
    elseif strcmp(Setting.choose_relax,'IRmax') % 指的是允许社会福利不够，尝试最大化IR(似乎没有意义)
        Cons = [Cons, Var.RelaxBB <= Setting.set_BB, Var.RelaxIC <= Setting.set_IC];
        obj = Var.RelaxIR; %最小化松弛量，相当于最大化盈余量
        Cons = [Cons, Var.RelaxIR >= -Setting.bigM];
        
    elseif strcmp(Setting.choose_relax,'ICmax') % 指的是允许社会福利不够，尝试最大化IC(似乎没有意义)
        Cons = [Cons, Var.RelaxBB <= Setting.set_BB, Var.RelaxIR <= Setting.set_IR];
        obj = Var.RelaxIC; %最小化松弛量，相当于最大化盈余量
        Cons = [Cons, Var.RelaxIC >= -Setting.bigM];
    end
    
    

    
    % 为什么不把这个约束写在basic problem里？因为Para.SW_max会变化
%     if isfield(Setting,'new_version') && Setting.new_version == 1 && Para.SW_max >0
%         Cons = [Cons, Var.SW == Para.SW_max * (1-Var.RelaxSW)];
%         Cons = [Cons, Var.BB == -Var.RelaxBB * Para.SW_max];
%     else
    Cons = [Cons, Para.SW_max - Var.SW - Setting.SWtolerance <= Var.RelaxSW];
%     end
%     Cons = [Cons, sum(sum(Para.Point .* Para.Prob .* Var.qip)) >= Para.SW_max - Var.RelaxSW - Setting.SWtolerance]; %社会福利. 只有这个社会福利约束最难达到

end
