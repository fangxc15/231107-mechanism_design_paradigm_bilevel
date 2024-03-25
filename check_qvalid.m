function [Resultout,solution] = check_qvalid(Result,Para,Num,Setting)

        %通过Qip来计算qip
        %根据场景集来编号
        Num.S = (Num.P+1) ^ Num.I;
        Var.Qipp = sdpvar(Num.I,Num.S,'full');
        Var.qip = sdpvar(Num.I,Num.P+1,'full');
        Var.qipadd = sdpvar(Num.I,Num.P+1,'full');
        Var.qipminus = sdpvar(Num.I,Num.P+1,'full');
        Var.qiporigin = sdpvar(Num.I, Num.P+1,'full');
        Cons = [];
        
        % 本身的容量上下限约束
        for i = 1:Num.I
            Cons = [Cons, Var.Qipp(i,:) >= Para.qmin(i)];
            Cons = [Cons, Var.Qipp(i,:) <= Para.qmax(i)];
        end
        
        % 本身的功率平衡约束
        Cons = [Cons, sum(Var.Qipp,1) == 0];
        Cons_num.Qbalance = length(Cons);
        
        % 现在可以再写一个拓扑约束. 拓扑约束必须要在Setting.topology=1的时候才成立
        if isfield(Setting,'topology')
            if Setting.topology == 1
            % 写出一个GSDF. GSDF是什么, Num.Branch*Num.I, Num.I*Num.S
                Cons = [Cons, Para.GSDF * Var.Qipp <= repmat(Para.BranchMax,1,Num.S)];
                Cons = [Cons, Para.GSDF * Var.Qipp >= repmat(-Para.BranchMax,1,Num.S)];
            end
        end
        
        
        scenario_prob = ones(1,Num.S);
        for i = 1:Num.I
            scenario_set(i,:) =  floor(mod(0:Num.S-1,(Num.P+1)^(Num.I- i +1))/((Num.P+1)^(Num.I-i)));
            scenario_prob =  scenario_prob .* Para.Prob(i,scenario_set(i,:)+1);
        end
        
        
        % 在qip和Qipp之间建立联系, 这个问题几乎肯定是infeasible
        for i = 1:Num.I
            for p = 1:Num.P+1
                Cons = [Cons, Var.qip(i,p) == sum(Var.Qipp(i,find(scenario_set(i,:)==p-1)) .* scenario_prob(find(scenario_set(i,:)==p-1)) / Para.Prob(i,p))];
                if abs(Para.Prob(i,p) - sum(scenario_prob(find(scenario_set(i,:)==p-1)))) > 1e-10
                    fprintf('fail %d %d',i,p)
                    disp(i)
                end 
                Cons_num.Qqlink(i,p) = length(Cons);
            end 
        end      
        
       for i = 1:Num.I
            % 容量上下限约束
            Cons = [Cons, Var.qip(i,:) >= Para.qmin(i) ];
            Cons = [Cons, Var.qip(i,:) <= Para.qmax(i) ];
       end 
       Cons = [Cons, Var.qip == Var.qiporigin + Var.qipadd - Var.qipminus];
       Cons = [Cons, Var.qipadd >=0, Var.qipminus >=0];
       
       Cons = [Cons, Var.qiporigin >= Result.qip];
       Dual.qip1 = length(Cons);
       Cons = [Cons, -Var.qiporigin >= -Result.qip];
       Dual.qip2 = length(Cons);

       obj = sum(sum(Var.qipadd + Var.qipminus));     

       ops = sdpsettings('solver','gurobi');
       solution = optimize(Cons,obj,ops);
       
       Resultout.obj = value(obj);
       % 计算总共的成交物品期望值
       Resultout.qip = value(Var.qip);
       Resultout.qipadd = value(Var.qipadd);
       Resultout.qipminus = value(Var.qipminus);
       Resultout.dualqip = dual(Cons(Dual.qip1)) - dual(Cons(Dual.qip2));
       Resultout.dualqip = reshape(Resultout.dualqip,Num.I,Num.P+1);
end

% 求解得到的这个没有意义,能施加约束吗?感觉不能
% 如果去求解它的对偶乘子,意义就是当qip扰动的时候可以减少多少的Qipp扰动。把这项放到主问题的决策里，也没有意义！我增加qip能减少多少Qipp的不平衡量，没意义
% 这个时候如果去添加割呢？ 添加割，是根据第二阶段的问题，去施加对第一阶段决策量的约束.

% 第一阶段的问题进行决策, 然后第二阶段的优化问题
% min 0
% 若干约束, AQ + Z= q, BQ<=Q_limit, CQ<=L_limit
% 原问题是否有解,其实真的不一定,
% 如果原问题有解,对偶问题的最优解是0，如果原问题无解,那对偶问题也可能无解或者是无界解.
% 这个可行问题有解(它肯定不是无界的),那它的对偶问题也一定要有解,且它的对偶问题的最优解一定要满足某个特征.
% 及w * [q, Q_limit, L_limit] >0, [A' B' C']  * w = 0(前面的0), 这是Benders约束.
% 具体如何添加一个合理的割呢? 因为根本就是无解, 哪来的w? 是否可以加入一个松弛项?

% 算法：
% 	求解主问题，得到第一阶段的决策变量y
% 	求解验证可行性子问题，若无解（不能让第二阶段子问题无解啊，不然的话哪来的u呢？除非求解对偶问题）
% 	可行性子问题是啥？把备用、功率两个地方的约束松弛一下。如果说要是解出来目标函数大于0，那说明有问题，也可以得到一个u。那说明要求这个u乘右手项后，小于0。作为一个割添加到主问题。
% 	再回到主问题，重新求解。
% 	最优性子问题
% 	如果求解最优性子问题，击穿了主问题里最优的 值，那么就新添加一条割。这里添加的割都是按照子问题对偶变量和子问题最优值一起写出的。
% 	击穿了，重新求解。


       


 
