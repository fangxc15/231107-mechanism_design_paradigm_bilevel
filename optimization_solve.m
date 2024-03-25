function [Result,solution] = optimization_solve(Num,Para,Setting, Var, Cons,obj)
    % 增加了通过不断添加割来简化单个问题，但最终完成迭代求解的方式
    % 但这里添加割有问题，有问题在于这个迭代次数太太太太多了. 那如何在一开始就预提交一些割? 这个是一个方法
    
    
    % 还有一个办法就是被打断的思路
    if Setting.direct_solve == 1
       
        ops = sdpsettings('bilevel.outersolver','bmibnb');
        % solvebilevel(CO,OO-OO^2,CI,OI,[y1 y2 y3],ops)
        Varset = [Var.QX];
        for c = 1:Num.C
            Varset = [Varset Var.scenario(c).QX Var.scenario(c).QXmean];
        end 
        solution = solvebilevel(Cons,obj,Cons_inner,-obj_inner,Varset,ops);
        Result.obj_inner = value(obj_inner); 
        Result.obj_base = value(obj_base);
        Result.obj_scenario = zeros(Num.C, Num.P);
        for c = 1:Num.C
            for p = 1:Num.P
                Result.obj_scenario(c,p) = value(obj_scenario(c,p).value);
            end
        end   
    else
        
        ops = sdpsettings('solver','gurobi');
        ops.verbose = 2;
%         ops = sdpsettings('solver','GUROBI','warning',0,'showprogress',1 ,'usex0',1) ; 
        ops.gurobi.MIPFocus = 0;      %默认值0,试图在最优值和可行解之间取得平衡。1：以可行解为目标；2：以最优解为目标；3：最优边界为目标。
        ops.gurobi.MIPGap = 0.001;   %默认值0.0001,从某Gap就停止，相对值
        ops.gurobi.MIRCuts = 2;       %Aggressive (2), Conservative (1), Automatic (-1), or None (0)
        ops.gurobi.Cuts = 3;          %0 to shut off cuts, 1 for moderate cut generation, 2 for aggressive cut generation, and 3 for very aggressive cut generation
        ops.gurobi.ImproveStartTime = 300;
        ops.gurobi.ImproveStartGap = 0.3;
        ops.gurobi.TimeLimit = 1200;
        ops.gurobi.Heuristics = 0.3;
        ops.cplex.MIPGap = 0.001; 
        
%         ops.gurobi.ResultFile = 'myfile.mps';
        ops.savesolveroutput = 1;
        ops.saveyalmipmodel = 1;
        ops.gurobi.TuneCriterion = 2;% TuneCriterion 调优准则。= 2 表示为目标函数值，即更专注于寻找好的可行解. 这是启发式算法的设计
%         ops.gurobi.Cutoff = 2.2924;%表示割集的上界要设定一下
        ops.gurobi.StartNodeLimit = 2000;% 如果传入初始解不完全，允许探索几个节点或许初始解，默认值是500
%         ops.gurobi.StartNumber %  从第几个初始解开始算
%         ops.gurobi.NumStart %有几个初始解
        
        
        if isfield(Setting,'Iter_init') && isfield(Setting,'Result_init')
            if Setting.Iter_init == 1
                Result_temp = Setting.Result_init; 
                if ~isempty(find(Setting.hard_assign==0))
%                 try
                                  
                    assign(Var.QX, Result_temp.QX);
                    for c = 1:Num.C
                        if Para.Cstrategic(c) == 1
                            assign(Var.scenario(c).QX, Result_temp.scenario(c).QX);
                            assign(Var.scenario(c).QXmean, Result_temp.scenario(c).QXmean);
                        end
                    end 

                    
                    for c = 1:Num.C
                        if Para.Cstrategic(c) == 1
                            assign(Var.scenario(c).dual_mumin, Result_temp.scenario(c).dual_mumin);
                            assign(Var.scenario(c).dual_mumax, Result_temp.scenario(c).dual_mumax);
                            assign(Var.scenario(c).dual_mumin_bin, Result_temp.scenario(c).dual_mumin_bin);
                            assign(Var.scenario(c).dual_mumax_bin, Result_temp.scenario(c).dual_mumax_bin);
                            assign(Var.scenario(c).dualprice, Result_temp.scenario(c).dualprice);
                            if Setting.topology == 1
                                assign(Var.scenario(c).dual_rhomax, Result_temp.scenario(c).dual_rhomax);
                                assign(Var.scenario(c).dual_rhomin, Result_temp.scenario(c).dual_rhomin);
                                assign(Var.scenario(c).dual_rhomax_bin, Result_temp.scenario(c).dual_rhomax_bin);
                                assign(Var.scenario(c).dual_rhomin_bin, Result_temp.scenario(c).dual_rhomin_bin);
                                assign(Var.scenario(c).Pl, Result_temp.scenario(c).Pl);
                                assign(Var.scenario(c).LMP, Result_temp.scenario(c).LMP);
                            end
                        end
                    end
                    assign(Var.dual_mumin, Result_temp.dual_mumin);
                    assign(Var.dual_mumax, Result_temp.dual_mumax);
                    assign(Var.dual_mumin_bin, Result_temp.dual_mumin_bin);
                    assign(Var.dual_mumax_bin, Result_temp.dual_mumax_bin);
                    assign(Var.dualprice, Result_temp.dualprice);
                    if Setting.topology == 1
                        assign(Var.dual_rhomax, Result_temp.dual_rhomax);
                        assign(Var.dual_rhomin, Result_temp.dual_rhomin);
                        assign(Var.dual_rhomax_bin, Result_temp.dual_rhomax_bin);
                        assign(Var.dual_rhomin_bin, Result_temp.dual_rhomin_bin);
                        assign(Var.Pl, Result_temp.Pl);
                        assign(Var.LMP, Result_temp.LMP);
                    end

                    assign(Var.spread, Result_temp.spread(:));
                    if Setting.unified_spread == 1    
                        assign(Var.unified_spread, Result_temp.unified_spread);
                    end


                    assign(Var.SW, Result_temp.SW);
                    assign(Var.Gcost, Result_temp.Gcost); 
                    assign(Var.Dutil, Result_temp.Dutil); 

                    if Setting.riskSW == 1
                        assign(Var.tempSW, Result_temp.tempSW);
                        assign(Var.baseSW, Result_temp.baseSW);
                        assign(Var.tempGcost, Result_temp.tempGcost);
                        assign(Var.baseGcost, Result_temp.baseGcost);
                        assign(Var.tempDutil, Result_temp.tempDutil);
                        assign(Var.baseDutil, Result_temp.baseDutil);
                    end
                    
                    
    
                    
                    assign(Var.MX, Result_temp.MX);
                    assign(Var.RX, Result_temp.RX);
                    assign(Var.BB, Result_temp.BB);
                    if Setting.topology == 1
                        assign(Var.congestion_revenue, Result_temp.congestion_revenue);
                    end
                    
                    assign(Var.RelaxBB, Result_temp.RelaxBB);
                    assign(Var.RelaxSW, Result_temp.RelaxSW);
                    assign(Var.RelaxIR, Result_temp.RelaxIR);
                    assign(Var.RelaxIC, Result_temp.RelaxIC);
                    
                
%                     for i = 1:length(Cons)
%                         primal_res(i) = check(Cons(i));
%                     end
%                     primal_res = primal_res(:);
%                     primal_res(abs(primal_res)<1e-10) = 0;

%                    for c = 1:Num.C 
%                        if Para.Cstrategic(c)==1
%                             temp = value(Var.scenario(c).QXmean(:,1) - (Var.scenario(c).QX(:,1) + Var.QX)/2);
%                             if max(max(temp))>0
%                                 fprintf([num2str(c),num2str(1)]);
%                                 temp
%                             end
%                             for p = 2:Num.P
%                                 temp = value(Var.scenario(c).QXmean(:,p) - (Var.scenario(c).QX(:,p) + Var.scenario(c).QX(:,p-1))/2);
%                                 if max(max(temp))>0
% 
%                                     fprintf([num2str(c),num2str(p)]);
% 
%                                 end
%                             end
%                         end
%                    end 
    


                    
                end
                if ~isempty(find(Setting.hard_assign==1))
                    
                    Cons = [Cons, Var.QX == Result_temp.QX];
                    for c =1:Num.C
                        if Para.Cstrategic(c) == 1
                            Cons = [Cons, Var.scenario(c).QX == Result_temp.scenario(c).QX];
                            Cons = [Cons, Var.scenario(c).QXmean == Result_temp.scenario(c).QXmean];
                        end
                    end 

                    
                    for c = 1:Num.C
                        if Para.Cstrategic(c) == 1
                            Cons = [Cons, Var.scenario(c).dual_mumin == Result_temp.scenario(c).dual_mumin];
                            Cons = [Cons, Var.scenario(c).dual_mumax == Result_temp.scenario(c).dual_mumax];
                            Cons = [Cons, Var.scenario(c).dualprice == Result_temp.scenario(c).dualprice];
                            Cons = [Cons, Var.scenario(c).dual_mumin_bin == Result_temp.scenario(c).dual_mumin_bin];
                            Cons = [Cons, Var.scenario(c).dual_mumax_bin == Result_temp.scenario(c).dual_mumax_bin];
                             
                            if Setting.topology == 1
                                Cons = [Cons, Var.scenario(c).dual_rhomax == Result_temp.scenario(c).dual_rhomax];
                                Cons = [Cons, Var.scenario(c).dual_rhomin == Result_temp.scenario(c).dual_rhomin];
                                Cons = [Cons, Var.scenario(c).dual_rhomax_bin == Result_temp.scenario(c).dual_rhomax_bin];
                                Cons = [Cons, Var.scenario(c).dual_rhomin_bin == Result_temp.scenario(c).dual_rhomin_bin];
                            end
                        end
                    end
                    Cons = [Cons, Var.dual_mumin == Result_temp.dual_mumin];
                    Cons = [Cons, Var.dual_mumax == Result_temp.dual_mumax];
                    Cons = [Cons, Var.dualprice == Result_temp.dualprice];
                    
                    Cons = [Cons, Var.dual_mumin_bin == Result_temp.dual_mumin_bin];
                    Cons = [Cons, Var.dual_mumax_bin == Result_temp.dual_mumax_bin];
                    
    
                    if Setting.topology == 1
                        Cons = [Cons, Var.dual_rhomax == Result_temp.dual_rhomax];
                        Cons = [Cons, Var.dual_rhomin == Result_temp.dual_rhomin];
                        Cons = [Cons, Var.dual_rhomax_bin == Result_temp.dual_rhomax_bin];
                        Cons = [Cons, Var.dual_rhomin_bin == Result_temp.dual_rhomin_bin];
                    end
                    
                    
%                     Cons = [Cons, Var.spread == Result_temp.spread(:)];
%                     
%                     if Setting.unified_spread == 1
%                         Cons = [Cons, Var.unified_spread == Result_temp.unified_spread];
%                     end
%                     end
                    Cons = [Cons, Var.SW == Result_temp.SW];
                    Cons = [Cons, Var.Gcost == Result_temp.Gcost];
                    Cons = [Cons, Var.Dutil == Result_temp.Dutil];
 

                    if Setting.riskSW == 1
                        Cons = [Cons, Var.tempSW == Result_temp.tempSW];
                        Cons = [Cons, Var.baseSW == Result_temp.baseSW];
                        Cons = [Cons, Var.tempGcost == Result_temp.tempGcost];
                        Cons = [Cons, Var.baseGcost == Result_temp.baseGcost];
                        Cons = [Cons, Var.tempDutil == Result_temp.tempDutil];
                        Cons = [Cons, Var.baseDutil == Result_temp.baseDutil];
                    end
                    
                    
    
                    Cons = [Cons, Var.MX == Result_temp.MX];
                    Cons = [Cons, Var.RX == Result_temp.RX];
                    Cons = [Cons, Var.BB == Result_temp.BB];
                    
                    
                    if Setting.topology == 1
                        Cons = [Cons, Var.congestion_revenue == Result_temp.congestion_revenue];
                    end
                    
                    Cons = [Cons, Var.RelaxBB == Result_temp.RelaxBB];
                    Cons = [Cons, Var.RelaxSW == Result_temp.RelaxSW];
                    Cons = [Cons, Var.RelaxIR == Result_temp.RelaxIR];
                    Cons = [Cons, Var.RelaxIC == Result_temp.RelaxIC];
     
                end
                
%                 ops = sdpsettings('solver','gurobi','usex0',1);
                ops.usex0 = 1;
                    
%                 catch
%                     fprintf(['There is no initial values for  market\n','\n']);
%                 end
                
                
            else
%                 ops = sdpsettings('solver','gurobi','usex0',0);
                ops.usex0 = 0;
            end
                
        end

        solution = optimize(Cons,obj,ops);
        
        Result.obj = value(obj);
        
        for c = 1:Num.C
            if Para.Cstrategic(c) == 1
                Result.scenario(c).dual_mumin = value(Var.scenario(c).dual_mumin);
                Result.scenario(c).dual_mumax = value(Var.scenario(c).dual_mumax);
                Result.scenario(c).dualprice = value(Var.scenario(c).dualprice); 
                Result.scenario(c).dual_mumin_bin = value(Var.scenario(c).dual_mumin_bin);
                Result.scenario(c).dual_mumax_bin = value(Var.scenario(c).dual_mumax_bin);
            end
        end
        Result.dual_mumin = value(Var.dual_mumin);
        Result.dual_mumax = value(Var.dual_mumax);
        Result.dual_mumin_bin = value(Var.dual_mumin_bin);
        Result.dual_mumax_bin = value(Var.dual_mumax_bin);
        Result.dualprice = value(Var.dualprice);
        if Setting.topology == 1
            Result.dual_rhomin = value(Var.dual_rhomin); 
            Result.dual_rhomax = value(Var.dual_rhomax);
            Result.dual_rhomin_bin = value(Var.dual_rhomin_bin); 
            Result.dual_rhomax_bin = value(Var.dual_rhomax_bin);
            
            
    %         Result.Pl =  -Para.GSDF * Result.QX;
            Result.Pl = value(Var.Pl);
            Result.LMP = value(Var.LMP);
%             Result.LMP = Result.dualprice + Para.GSDF' * (Result.dual_rhomin - Result.dual_rhomax);
            %LMP推导方式，在min问题中
            %sum QX + D <=0
            %- Para.GSDF * QX - Para.GSDF * D<=Para.BranchMax
            %这个LMP原来写的方式是相反的
            for c = 1:Num.C
                if Para.Cstrategic(c) == 1
                    Result.scenario(c).dual_rhomin = value(Var.scenario(c).dual_rhomin);
                    Result.scenario(c).dual_rhomax = value(Var.scenario(c).dual_rhomax);
                    Result.scenario(c).dual_rhomin_bin = value(Var.scenario(c).dual_rhomin_bin);
                    Result.scenario(c).dual_rhomax_bin = value(Var.scenario(c).dual_rhomax_bin);
                end
            end
        else
            Result.LMP = repmat(Result.dualprice,Num.I,1);
        end
%         Result.LMPsettle = Result.dualprice * Result.QX; %我真的是要按照LMP settle吗？我觉得不是的，我应该分别按照买方,卖方的LMP settle.     
    end
    
%     Result.obj = value(obj);
    
    Result.QX = value(Var.QX);
    for c = 1:Num.C 
        if Para.Cstrategic(c) == 1
            Result.scenario(c).QX = value(Var.scenario(c).QX);
            Result.scenario(c).QXmean = value(Var.scenario(c).QXmean);
            if Setting.topology == 1
                Result.scenario(c).Pl = value(Var.scenario(c).Pl);
                Result.scenario(c).LMP = value(Var.scenario(c).LMP);
            end
        end
    end 
    
    if Setting.unified_spread == 1    
        Result.unified_spread = value(Var.unified_spread);
    end
    Result.spread = value(Var.spread);
    
    Result.SW = value(Var.SW);
    Result.Gcost = value(Var.Gcost);
    Result.Dutil = value(Var.Dutil);

    if Setting.riskSW == 1
        Result.tempSW = value(Var.tempSW);
        Result.tempGcost = value(Var.tempGcost);
        Result.tempDutil = value(Var.tempDutil);
        Result.baseSW = value(Var.baseSW);
        Result.baseGcost = value(Var.baseGcost);
        Result.baseDutil = value(Var.baseDutil);
    end

    Result.MX = value(Var.MX);
    Result.RX = value(Var.RX);
    Result.BB = value(Var.BB);


    
    Result.RelaxBB = value(Var.RelaxBB);
    Result.RelaxSW = value(Var.RelaxSW);
    Result.RelaxIR = value(Var.RelaxIR);
    Result.RelaxIC = value(Var.RelaxIC);
    
    if Setting.riskSW == 1
        Result.qtotal =  sum(abs(Result.QX)) * Para.prob(1,1);
        for c = 1:Num.C
            if Para.Cstrategic(c) == 1
                for p = 1:Num.P
                    Result.qtotal = Result.qtotal + Para.prob(c,p+1) * sum(abs(Result.scenario(c).QX(:,p)));
                end
            end
        end  
    else
        Result.qtotal = sum(abs(Result.QX));
    end
    Result.welfare = Result.SW;
    Result.surplus = Result.BB;
    
    if Setting.direct_solve == 0
        Result.LMPI = (Result.LMP + Para.signI' .* Result.spread);
        Result.LMPsettle = Result.LMPI .* Result.QX;
        Result.LMPsettleC = zeros(Num.C,1);
        for c = 1:Num.C
            Result.LMPsettleC(c) = sum(Result.LMPsettle(Para.company_set(c).no));
        end
        Result.Settle_diff = Result.MX - Result.LMPsettleC;
    end
    
    
    if Setting.topology == 1 && Setting.direct_solve == 0
        Result.congestion_revenue = value(Var.congestion_revenue);
        Result.LMPbalance = sum(Result.LMPsettleC) - Result.congestion_revenue;
    end
    
    Result.kIC = value(Var.kIC);
    Result.kIR = value(Var.kIR);
    Result.RelaxIR_C = value(Var.RelaxIR_C);
    if Setting.riskSW == 0
        Result.clearing_obj = (Para.xpara - Para.sign.* Result.spread(:)') * Result.QX;
    else     
        Result.tempobj = zeros(Num.C,Num.P);
        for c = 1:Num.C
            if Para.Cstrategic(c) == 1
                for p = 1:Num.P
                    temppara = Para.xpara;
                    tempno = Para.company_set(c).no;
                    temppara(tempno) = Para.Point(tempno, p+1);
                    Result.tempobj(c,p) = (temppara  - Para.signI.* Result.spread(:)') * Result.scenario(c).QX(:,p);
                end
            end
        end
        Result.baseobj = (Para.xpara  - Para.signI.* Result.spread(:)') * Result.QX;        
        Result.clearing_obj = Result.baseobj + sum(sum(Result.tempobj));      
%         Cons = [Cons, Var.tempobj(c,p) == (temppara  - Para.signI.* input_spread) * Var.scenario(c).QX(:,p)];
%         Cons = [Cons, Var.baseobj == (Para.xpara  - Para.signI.* input_spread) * Var.QX];
%         Cons = [Cons, Var.obj == Var.baseobj + sum(sum(Var.tempobj))];
    end
    
    Result.QXPara = zeros(Num.C,1);
    for c = 1:Num.C %这个可以补充进来, 但这里指的是基础状态下的生产成本
        tempno = Para.company_set(c).no;
        Result.QXPara(c) = Para.xpara(tempno) * Result.QX(tempno);       
    end
    Result.producer_welfare = sum(Result.RX(Para.Gset_company));
    Result.consumer_welfare = sum(Result.RX(Para.Dset_company));
    if Setting.direct_solve == 0
        Result.LMPsettleC_RX = Result.QXPara - Result.LMPsettleC;
    end
    




%     sum(Setting.Result_primal.LMPsettleC) + Setting.Result_primal.congestion_revenue
    
%     value(Para.signI' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin + ...
%                                                 Para.GSDF' * (Var.dual_rhomax - Var.dual_rhomin) + Var.dualprice)
%     value(Para.signI' .* Var.spread)
%     Para.xpara'
%     value(Var.dual_mumax)
%     value(Var.dual_mumin)
%     value(Var.dualprice)
%     value(Para.GSDF' * (Var.dual_rhomax - Var.dual_rhomin))
                                            




%     Result.obj = value(obj);
%     % 计算总共的成交物品期望值
%     Result.qtotal = sum(sum(Para.Prob(Para.Dset,:) .* value(Var.qip(Para.Dset,:))));
%     Result.Gcost = sum(sum(Para.Point(Para.Gset,:) .* Para.Prob(Para.Gset,:) .* value(Var.qip(Para.Gset,:))));
%     Result.Duti = sum(sum(Para.Point(Para.Dset,:) .* Para.Prob(Para.Dset,:) .* value(Var.qip(Para.Dset,:))));
%     Result.Grevenue = sum(sum(Para.Prob(Para.Gset,:) .* value(Var.mip(Para.Gset,:))));
%     Result.Dbill = sum(sum(Para.Prob(Para.Dset,:) .* value(Var.mip(Para.Dset,:))));
%         % 确认物品是否能够严格Balance
%     Result.qbalance = value(sum(sum(Para.Prob .* Var.qip)));
%     Result.welfare = value(sum(sum(Para.Point .* Para.Prob .* Var.qip)));
%     Result.surplus = value(sum(sum(Para.Prob .* Var.mip)));
    
%     Result.RelaxBB = value(Var.RelaxBB);
%     Result.RelaxSW = value(Var.RelaxSW);
%     Result.RelaxIR = value(Var.RelaxIR); %这两个还有可以修改的地方。松弛IC后，IR也不是一定可以满足的哦！
%     Result.RelaxIC = value(Var.RelaxIC); 
%     fprintf('物品平衡')
%     disp(Result.qbalance)
%     fprintf('社会福利')
%     disp(Result.welfare)
%     fprintf('收支盈余量')
%     disp(Result.surplus)
%     fprintf('个体理性松弛量')
%     disp(Result.RelaxIR)
%     fprintf('激励相容松弛量')
%     disp(Result.RelaxIC)
% 
% 
% 
% 
%     % 我这里既然写出了IR(贝叶斯意义上的)和IC,那我也可以根据定义重新写出一个IR.
%     Result.qip = value(Var.qip);
%     Result.mip = value(Var.mip);
%     Result.rip = value(Var.rip);
%     Result.qip_mean = value(Var.qip_mean);
%     if isfield(Setting,'topology') && Setting.Qipp == 1
%         if Setting.topology == 1
%             Result.Lipp = Para.GSDF * value(Var.Qipp);
%         end
%     end
    
end

