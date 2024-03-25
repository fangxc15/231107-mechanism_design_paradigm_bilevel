function [Result,solution] = optimization_solve_origin(Num,Para,Setting, Var, Cons,obj)
    % 原本里没有考虑当Qipp=0的时候的求解方式
    ops = sdpsettings('solver','gurobi');


    solution = optimize(Cons,obj,ops);


    Result.obj = value(obj);
    % 计算总共的成交物品期望值
    Result.qtotal = sum(sum(Para.Prob(Para.Dset,:) .* value(Var.qip(Para.Dset,:))));
    Result.Gcost = sum(sum(Para.Point(Para.Gset,:) .* Para.Prob(Para.Gset,:) .* value(Var.qip(Para.Gset,:))));
    Result.Duti = sum(sum(Para.Point(Para.Dset,:) .* Para.Prob(Para.Dset,:) .* value(Var.qip(Para.Dset,:))));
    Result.Grevenue = sum(sum(Para.Prob(Para.Gset,:) .* value(Var.mip(Para.Gset,:))));
    Result.Dbill = sum(sum(Para.Prob(Para.Dset,:) .* value(Var.mip(Para.Dset,:))));
    fprintf('成交总量')
    disp(Result.qtotal)
    fprintf('发电成本')
    disp(Result.Gcost)
    fprintf('用户效用')
    disp(Result.Duti)
    fprintf('发电收入')
    disp(Result.Grevenue)
    fprintf('用电支出')
    disp(Result.Dbill)


    % 确认物品是否能够严格Balance
    Result.qbalance = value(sum(sum(Para.Prob .* Var.qip)));
    Result.welfare = value(sum(sum(Para.Point .* Para.Prob .* Var.qip)));
    Result.surplus = value(sum(sum(Para.Prob .* Var.mip)));
    
    Result.RelaxBB = value(Var.RelaxBB);
    Result.RelaxSW = value(Var.RelaxSW);
    Result.RelaxIR = value(Var.RelaxIR); %这两个还有可以修改的地方。松弛IC后，IR也不是一定可以满足的哦！
    Result.RelaxIC = value(Var.RelaxIC); 
    fprintf('物品平衡')
    disp(Result.qbalance)
    fprintf('社会福利')
    disp(Result.welfare)
    fprintf('收支盈余量')
    disp(Result.surplus)
    fprintf('个体理性松弛量')
    disp(Result.RelaxIR)
    fprintf('激励相容松弛量')
    disp(Result.RelaxIC)




    % 我这里既然写出了IR(贝叶斯意义上的)和IC,那我也可以根据定义重新写出一个IR.
    Result.qip = value(Var.qip);
    Result.mip = value(Var.mip);
    Result.rip = value(Var.rip);
    Result.qip_mean = value(Var.qip_mean);
    if isfield(Setting,'topology') && Setting.Qipp == 1
        if Setting.topology == 1
            Result.Lipp = Para.GSDF * value(Var.Qipp);
        end
    end
end