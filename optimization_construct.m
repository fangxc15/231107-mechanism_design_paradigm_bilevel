function  [Var, Cons, Cons_num] = optimization_construct(Num,Para,Setting)
    
    % 在这里，我甚至无需去设计Q(X)单调递增的约束, 因为它肯定是单调递增的.
    % 要注意互补松弛条件写的时候Para.bigM的数值问题
    Cons_num = 0;
    Var.spread = sdpvar(Num.I,1); %这个spread可以是对所有主体的
    Var.kIC = sdpvar(1,1); 
    Var.kIR = sdpvar(1,1);
    Var.RelaxIR_C = sdpvar(Num.C,1); % 分主体的relax_IR
    
    
    Var.MX = sdpvar(Num.C,1); 
    Var.RX = sdpvar(Num.C,1); 

    Var.RelaxBB = sdpvar(1,1);
    Var.RelaxSW = sdpvar(1,1);
    Var.RelaxIC = sdpvar(1,1);
    Var.RelaxIR = sdpvar(1,1);
    Var.BB = sdpvar(1,1);
    Var.SW = sdpvar(1,1);
    Var.Gcost = sdpvar(1,1);
    Var.Dutil = sdpvar(1,1);
%     Var.Grevenue = sdpvar(1,1);
%     Var.Drevenue = sdpvar(1,1);

    % Variables for lower level
    Cons = [];
    Var.QX = sdpvar(Num.I,1);
    for c = 1:Num.C 
        if Para.Cstrategic(c) == 1
            Var.scenario(c).QX = sdpvar(Num.I,Num.P,'full');
            Var.scenario(c).QXmean = sdpvar(Num.I,Num.P,'full');
            % 想要表达的约束是，company自己随着自己的参数变化都
            tempno = Para.company_set(c).no;
            % 这个约束表示QX>QX_scenario. 想要表达的是随着参数X增加QX单调增.
            % 如果是卖方主体，从x一直到x_worse, interval_len < 0, 在增加的过程中，前者减后者也小于0
            Cons = [Cons, Para.Interval_len(tempno,1)' * (Var.QX(tempno) - Var.scenario(c).QX(tempno,1)) >= 0];
            for p = 2:Num.P
                Cons = [Cons, Para.Interval_len(tempno,p)' * (Var.scenario(c).QX(tempno,p-1) - Var.scenario(c).QX(tempno,p))>= 0];
                    %对于买方来讲，参数越好，中标量越大, QX越大
                    %对于卖方来讲，参数越好，中标量越大，QX越小
            end   
        end
    end 
    

    
    Cumulative_cons = 0;
    % Constraints on spread design
    Cons = [Cons, Var.spread >= 0];
    Cons = [Cons, Var.spread <= max(Para.xmax)/2];
    
    Cumulative_cons = Cumulative_cons +  Num.I;
    
    
    if Setting.unified_spread == 1    
        Var.unified_spread = sdpvar(1,1);
        Cons = [Cons, Var.spread == Var.unified_spread]; Cumulative_cons = Cumulative_cons +  Num.I;
        
    end 

    % Calculate Welfare
    if Setting.riskSW == 0
        Cons = [Cons, Para.xpara(Para.Gset) * Var.QX(Para.Gset) == Var.Gcost]; Cumulative_cons = Cumulative_cons +  1;
        Cons = [Cons, Para.xpara(Para.Dset) * Var.QX(Para.Dset) == Var.Dutil]; Cumulative_cons = Cumulative_cons +  1;
        Cons = [Cons, Var.Gcost + Var.Dutil == Var.SW];                        Cumulative_cons = Cumulative_cons +  1;
%         Cons = [Cons, Var.Grevenue == Var.Gcost + sum(Var.MX(Para.Gset))];
%         Cons = [Cons, Var.Drevenue == Var.Dutil + sum(Var.MX(Para.Dset))];
    else
        Var.tempSW = sdpvar(Num.C, Num.P,'full');
        Var.tempGcost = sdpvar(Num.C, Num.P,'full');
        Var.tempDutil = sdpvar(Num.C, Num.P,'full');
        Var.baseSW = sdpvar(1,1);
        Var.baseGcost = sdpvar(1,1);
        Var.baseDutil = sdpvar(1,1);

        for c = 1:Num.C
            if Para.Cstrategic(c) == 0
%                 Cons = [Cons, Var.scenario(c).QX == 0]; Cumulative_cons = Cumulative_cons +  Num.I;
                Cons = [Cons, Var.tempSW(c,:) == 0];
                Cons = [Cons, Var.tempGcost(c,:) == 0];
                Cons = [Cons, Var.tempDutil(c,:) == 0];
            else
                for p = 1:Num.P
                    temppara = Para.xpara;
                    tempno = Para.company_set(c).no;
                    temppara(tempno) = Para.Point(tempno, p+1);
                    Cons = [Cons, temppara * Var.scenario(c).QX(:,p) == Var.tempSW(c,p)]; Cumulative_cons = Cumulative_cons +  Num.I;
                    Cons = [Cons, temppara(Para.Gset) * Var.scenario(c).QX(Para.Gset,p) == Var.tempGcost(c,p)]; Cumulative_cons = Cumulative_cons +  1;
                    Cons = [Cons, temppara(Para.Dset) * Var.scenario(c).QX(Para.Dset,p) == Var.tempDutil(c,p)]; Cumulative_cons = Cumulative_cons +  1;
                end
            end
        end
        Cons = [Cons, Para.xpara(Para.Gset) * Var.QX(Para.Gset) == Var.baseGcost]; Cumulative_cons = Cumulative_cons +  1;
        Cons = [Cons, Para.xpara(Para.Dset) * Var.QX(Para.Dset) == Var.baseDutil]; Cumulative_cons = Cumulative_cons +  1;
        Cons = [Cons, Var.baseSW == Para.xpara * Var.QX];                          Cumulative_cons = Cumulative_cons +  1;
        
        Cons = [Cons, Var.SW/Para.prob(1,1) == Var.baseSW + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempSW))/Para.prob(1,1)];          Cumulative_cons = Cumulative_cons +  1;
        Cons = [Cons, Var.Gcost/Para.prob(1,1) == Var.baseGcost + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempGcost))/Para.prob(1,1)]; Cumulative_cons = Cumulative_cons +  1;
        Cons = [Cons, Var.Dutil/Para.prob(1,1) == Var.baseDutil + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempDutil))/Para.prob(1,1)]; Cumulative_cons = Cumulative_cons +  1;

%         Cons = [Cons, Var.SW == Para.prob(1,1) * Var.baseSW + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempSW))];          Cumulative_cons = Cumulative_cons +  1;
%         Cons = [Cons, Var.Gcost == Para.prob(1,1) * Var.baseGcost + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempGcost))]; Cumulative_cons = Cumulative_cons +  1;
%         Cons = [Cons, Var.Dutil == Para.prob(1,1) * Var.baseDutil + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempDutil))]; Cumulative_cons = Cumulative_cons +  1;
    end

    % Calculate QX mean
    for c = 1:Num.C 
        if Para.Cstrategic(c) == 1
            Cons = [Cons, Var.scenario(c).QXmean(:,1) == (Var.scenario(c).QX(:,1) + Var.QX)/2];                          Cumulative_cons = Cumulative_cons +  Num.I;
            for p = 2:Num.P
                Cons = [Cons, Var.scenario(c).QXmean(:,p) == (Var.scenario(c).QX(:,p) + Var.scenario(c).QX(:,p-1))/2];   Cumulative_cons = Cumulative_cons +  Num.I;
            end
        end
    end 

    % Tempobj
%     obj = Var.RelaxSW;
%     Cons = [Cons, Var.RelaxIC <=0, Var.RelaxIR <= 0,Var.RelaxBB <=0];

    %% Direct Design of Lower Level Model
    if Setting.direct_solve == 1
        % Lower Level Model Design
        obj_base = (Para.xpara  - Para.sign .* Var.spread')* Var.QX; % Clearing for one case
        obj_inner = obj_base;
        Cons_inner = [];
        Cons_inner = [Cons_inner, sum(Var.QX) == 0];
        Cons_inner = [Cons_inner, Var.QX <= Para.qmax', Var.QX >= Para.qmin'];
        if Setting.topology == 1
            % 写出一个GSDF. GSDF是什么, Num.Branch*Num.I, Num.I*Num.S
            Cons = [Cons, -Para.GSDF * Var.QX <= Para.BranchMax];
            Cons = [Cons, -Para.GSDF * Var.QX >= Para.BranchMax];
        end
        for c = 1:Num.C 
            for p = 1:Num.P
                temppara = Para.xpara;
                tempno = Para.company_set(c).no;
                temppara(tempno) = Para.Point(tempno, p+1);
                obj_scenario(c,p).value = (temppara - Para.sign .* Var.spread')* Var.scenario(c).QX(:,p); % Clearing for one case
                obj_inner = obj_inner + obj_scenario(c,p).value;
                Cons_inner = [Cons_inner, sum(Var.scenario(c).QX(:,p)) == 0];
                Cons_inner = [Cons_inner, Var.scenario(c).QX(:,p) <= Para.qmax', Var.scenario(c).QX(:,p) >= Para.qmin'];
            end 
        end
    elseif Setting.direct_solve == 0
        
        Var.dualprice = sdpvar(1,1);
        Var.dual_mumax = sdpvar(Num.I,1);
        Var.dual_mumin = sdpvar(Num.I,1);
        Var.dual_mumax_bin = binvar(Num.I,1);
        Var.dual_mumin_bin = binvar(Num.I,1);
        
        Cons_inner = [];
        Cons_inner = [Cons_inner, sum(Var.QX) == 0]; Cumulative_cons = Cumulative_cons +  Num.I;
        % for i = 1:Num.I
        Cons_inner = F_Complementarity_Yal(Cons_inner, Var.QX - Para.qmin', Var.dual_mumin, Var.dual_mumin_bin, Setting);
        Cumulative_cons = Cumulative_cons +  4*Num.I;
        Cons_inner = F_Complementarity_Yal(Cons_inner, Para.qmax' - Var.QX, Var.dual_mumax, Var.dual_mumax_bin, Setting);
        Cumulative_cons = Cumulative_cons +  4*Num.I;
        if Setting.topology == 1     
            Var.dual_rhomax = sdpvar(Num.Branch,1);
            Var.dual_rhomin = sdpvar(Num.Branch,1);
            Var.dual_rhomax_bin = binvar(Num.Branch,1);
            Var.dual_rhomin_bin = binvar(Num.Branch,1);
            
            Var.Pl = sdpvar(Num.Branch,1);
            Var.LMP = sdpvar(Num.I,1);

            Cons_inner = [Cons_inner, 100 * Para.GSDF * Var.QX == -100*Var.Pl];
            Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax - Var.Pl , Var.dual_rhomax, Var.dual_rhomax_bin, Setting);
%             Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax + Para.GSDF * Var.QX , Var.dual_rhomax, Var.dual_rhomax_bin, Setting);
            Cumulative_cons = Cumulative_cons +  4*Num.Branch;
            Cons_inner = F_Complementarity_Yal(Cons_inner, Var.Pl + Para.BranchMax,  Var.dual_rhomin, Var.dual_rhomin_bin, Setting);
%             Cons_inner = F_Complementarity_Yal(Cons_inner, - Para.GSDF * Var.QX + Para.BranchMax,  Var.dual_rhomin, Var.dual_rhomin_bin, Setting);
            Cumulative_cons = Cumulative_cons +  4*Num.Branch;
            
            
            
            Cons_inner = [Cons_inner,  100 * Var.LMP == 100 * Var.dualprice - 100 * Para.GSDF' * (Var.dual_rhomax - Var.dual_rhomin)];
            Cons_inner = [Cons_inner, Para.signI' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin + Var.LMP == 0];
%             Cons_inner = [Cons_inner, Para.signI' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin - ...
%                                                 Para.GSDF' * (Var.dual_rhomax - Var.dual_rhomin) + Var.dualprice == 0];
            Cumulative_cons = Cumulative_cons +  Num.I;
            % QX>0表示买家获得能量, max问题要转换为min问题, min -[Para +
            % spread(买家为-/卖家为+)], 所有的小于等于约束前都为正。这里的Para.Branch疑似要去反过来. 功率平衡约束，本质上是一个<=0的约束，需求5 供给-7，要小于等于0]
        else
            %目前暂时的不含有拓扑的出清,Lagrange条件
            
            Cons_inner = [Cons_inner, Para.signI' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin + Var.dualprice == 0];
            Cumulative_cons = Cumulative_cons +  Num.I;
        end
        
        
        
        for c = 1:Num.C
            if Para.Cstrategic(c) == 1
                Var.scenario(c).dualprice = sdpvar(1,Num.P,'full');
                Var.scenario(c).dual_mumax = sdpvar(Num.I,Num.P,'full');
                Var.scenario(c).dual_mumin = sdpvar(Num.I,Num.P,'full');
                Var.scenario(c).dual_mumax_bin = binvar(Num.I,Num.P,'full');
                Var.scenario(c).dual_mumin_bin = binvar(Num.I,Num.P,'full');

                if Setting.topology == 1
                    Var.scenario(c).dual_rhomax = sdpvar(Num.Branch,Num.P,'full');
                    Var.scenario(c).dual_rhomin = sdpvar(Num.Branch,Num.P,'full');
                    Var.scenario(c).dual_rhomax_bin = binvar(Num.Branch,Num.P,'full');
                    Var.scenario(c).dual_rhomin_bin = binvar(Num.Branch,Num.P,'full');
                    
                    Var.scenario(c).LMP = sdpvar(Num.I,Num.P,'full');
                    Var.scenario(c).Pl = sdpvar(Num.Branch,Num.P,'full');
                end
            
            
                Cons_inner = [Cons_inner, sum(Var.scenario(c).QX,1) == 0]; Cumulative_cons = Cumulative_cons +  Num.I;
                Cons_inner = F_Complementarity_Yal(Cons_inner, Var.scenario(c).QX - repmat(Para.qmin',1,Num.P), Var.scenario(c).dual_mumin, Var.scenario(c).dual_mumin_bin, Setting);
                Cumulative_cons = Cumulative_cons +  4 * Num.I;
                Cons_inner = F_Complementarity_Yal(Cons_inner, repmat(Para.qmax',1,Num.P) - Var.scenario(c).QX, Var.scenario(c).dual_mumax, Var.scenario(c).dual_mumax_bin, Setting);
                Cumulative_cons = Cumulative_cons +  4 * Num.I;
    %                 for p = 1:Num.P
    %                     Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax - Para.GSDF * Var.scenario(c).QX(:,p) ,...
    %                         Var.scenario(c).dual_rhomax(:,p), Var.scenario(c).dual_rhomax_bin(:,p), Setting);
    %                     Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax + Para.GSDF * Var.scenario(c).QX(:,p) ,...
    %                         Var.scenario(c).dual_rhomin(:,p), Var.scenario(c).dual_rhomin_bin(:,p), Setting);
    %             
    %                 Cons_inner = [Cons_inner, Para.signI' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin + ...
    %                                                 Para.GSDF' * (Var.dual_rhomax - Var.dual_rhomin) + Var.dualprice == 0];
    %                 
    %             end

                for p = 1:Num.P
                    temppara = Para.xpara;
                    tempno = Para.company_set(c).no;
                    temppara(tempno) = Para.Point(tempno, p+1);
                    if Setting.topology == 1
                        
                        Cons_inner = [Cons_inner, 100 * Para.GSDF * Var.scenario(c).QX(:,p) == - 100 * Var.scenario(c).Pl(:,p)];
                        Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax - Var.scenario(c).Pl(:,p) ,...
                            Var.scenario(c).dual_rhomax(:,p), Var.scenario(c).dual_rhomax_bin(:,p), Setting);
                        Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax + Var.scenario(c).Pl(:,p) ,...
                            Var.scenario(c).dual_rhomin(:,p), Var.scenario(c).dual_rhomin_bin(:,p), Setting);

%                         Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax + Para.GSDF * Var.scenario(c).QX(:,p) ,...
%                             Var.scenario(c).dual_rhomax(:,p), Var.scenario(c).dual_rhomax_bin(:,p), Setting);
%                         Cons_inner = F_Complementarity_Yal(Cons_inner, Para.BranchMax - Para.GSDF * Var.scenario(c).QX(:,p) ,...
%                             Var.scenario(c).dual_rhomin(:,p), Var.scenario(c).dual_rhomin_bin(:,p), Setting);
                        Cons_inner = [Cons_inner, 100 * Var.scenario(c).LMP(:,p) == - 100 * Para.GSDF' * (Var.scenario(c).dual_rhomax(:,p) - ...
                                            Var.scenario(c).dual_rhomin(:,p)) + 100 * repmat(Var.scenario(c).dualprice(1,p),Num.I,1)];
                        Cons_inner = [Cons_inner, Para.signI' .* Var.spread - temppara' + Var.scenario(c).dual_mumax(:,p) ... 
                                            - Var.scenario(c).dual_mumin(:,p) + Var.scenario(c).LMP(:,p) == 0];
                
            

%                         Cons_inner = [Cons_inner, Para.signI' .* Var.spread - temppara' + Var.scenario(c).dual_mumax(:,p) ... 
%                                 - Var.scenario(c).dual_mumin(:,p) - Para.GSDF' * (Var.scenario(c).dual_rhomax(:,p) - ...
%                                 Var.scenario(c).dual_rhomin(:,p)) + repmat(Var.scenario(c).dualprice(1,p),Num.I,1) == 0];
                    else
                        Cons_inner = [Cons_inner, Para.signI' .* Var.spread - temppara' + Var.scenario(c).dual_mumax(:,p) ... 
                                                    - Var.scenario(c).dual_mumin(:,p) + repmat(Var.scenario(c).dualprice(1,p),Num.I,1) == 0];
                    end
                end
            end
        end 
        Cons = [Cons, Cons_inner];
    end
    
    Para.qabsmax = zeros(1,Num.I);
    Para.qabsmax(Para.Dset) = Para.qmax(Para.Dset);
    Para.qabsmax(Para.Gset) = abs(Para.qmin(Para.Gset));
    % Calculate MX
    for c = 1:Num.C
        tempno = Para.company_set(c).no;
%         Para.qabsmaxC(c) = sum(Para.qabsmax(tempno));
        if Para.Cstrategic(c) == 1
            % 如果是策略性的主体
            if length(tempno) > 1 % a company with multiple units
                temp_len = sum(sum(Para.Interval_len(tempno,:) .^2,1).^0.5); %这里的temp_len是总共的len
                if isfield(Setting,'new_version') && Setting.new_version == 1
                    % 这里引入kIR和kIC分摊
                    
                    % 假设说是卖方，此时的interval_len <0 (从x - x_worse), 而此时的QXmean<0，相乘后大于0
                    % 假设说是买方，此时的interval_len >0 (从x - x_worse),
                    % 而此时的QXmean>0, 相乘后大于0.
                    Cons = [Cons, Var.RX(c,1) == sum(sum(Para.Interval_len(tempno,:) .* ...
                        Var.scenario(c).QXmean(tempno,:),1)) - Var.kIC * temp_len * ...
                        sum(Para.qabsmax(tempno)) - Var.kIR * sum(Para.qabsmax(tempno))];
                elseif  isfield(Setting,'new_version') && Setting.new_version == 2
                    % 这里是kIC分摊的时候还要Max(0)
                    temp_len = sum(Para.Interval_len(tempno,:) .^2,1).^0.5;
                     Cons = [Cons, Var.RX(c,1) == sum(max(0,sum(Para.Interval_len(tempno,:) .* ...
                        Var.scenario(c).QXmean(tempno,:),1) - Var.kIC * sum(Para.qabsmax(tempno)) ...
                        * temp_len))  - Var.kIR * sum(Para.qabsmax(tempno))];
                elseif  isfield(Setting,'new_version') && Setting.new_version == 11
                    % 这里是改变kIC定义的方式，由统一乘一个templen变为分段相乘后求和. sumproduct
                     Cons = [Cons, Var.RX(c,1) == sum(max(0,sum(Para.Interval_len(tempno,:) .* ...
                        Var.scenario(c).QXmean(tempno,:),1) - Var.kIC * abs(Para.qabsmax(tempno)) * abs(Para.Interval_len(tempno,:)) ...
                        ))  - Var.kIR * sum(Para.qabsmax(tempno))];
                    
                    
                else
                    % 最早直接按relaxIR进行分摊
                    Cons = [Cons, Var.RX(c,1) == sum(sum(Para.Interval_len(tempno,:) .* ...
                        Var.scenario(c).QXmean(tempno,:),1)) - Var.RelaxIC * temp_len - Var.RelaxIR];
                end
                Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Var.RX(c,1)];
            else % a single unit
                if isfield(Setting,'new_version') && Setting.new_version == 1
                    % 这里如果是买方主体，此时的Interval_len>0, sign>0, 此时想办法减去kIC部分
                    % 这里如果是卖方主体，此时的Interval_len<0, sign<0, 还是想办法减去kIC部分
                    Cons = [Cons, Var.RX(c,1) == sum(Para.Interval_len(tempno,:) .* ...
                        (Var.scenario(c).QXmean(tempno,:) - Para.sign(c) * Var.kIC * ...
                        sum(Para.qabsmax(tempno)) )) - Var.kIR * sum(Para.qabsmax(tempno))];
                elseif  isfield(Setting,'new_version') && (Setting.new_version == 2 || Setting.new_version == 11)    
                    Cons = [Cons, Var.RX(c,1) == sum(max(0,Para.Interval_len(tempno,:) .* ...
                        (Var.scenario(c).QXmean(tempno,:) - Para.sign(c) * Var.kIC * ...
                        sum(Para.qabsmax(tempno)) ))) - Var.kIR * sum(Para.qabsmax(tempno))];
                else
                     Cons = [Cons, Var.RX(c,1) == sum(Para.Interval_len(tempno,:) .* ...
                        (Var.scenario(c).QXmean(tempno,:) - Para.sign(c) * Var.RelaxIC)) - Var.RelaxIR];
                end
                Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Var.RX(c,1)];
    %             Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - sum(Para.Interval_len(tempno,:) .* ...
    %                 (Var.scenario(c).QXmean(tempno,:) - Para.sign(c) * Var.RelaxIC)) + Var.RelaxIR];
            end
        else % Ruguo不是策略性主体
            % 这里应该会有LMP和QX相乘，有非线性项
            % 已知LMP = Result.dualprice + Para.GSDF' * (Result.dual_rhomin - Result.dual_rhomax);
            % 又根据Lagrange条件已知
            % Para.signI' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin + ...
            %           Para.GSDF' * (Var.dual_rhomin - Var.dual_rhomax) + Var.dualprice == 0
            % LMP = - Para.signI' .* Var.spread + Para.xpara' - Var.dual_mumax + Var.dual_mumin
            % 此时LMP和QX相乘，关键难点在于Var.spread 和 QX相乘
            % 根据底层问题的强对偶定理
            % \sum (Para.xpara +- Var.spread) * Var.QX
            % 如果想要求某个Var.spread * Var.QX,其实是办不到的。
            
            % 当然，其实真正结算的价格是LMP + Para.signI' .* Var.spread, 需要将QX和这一项相乘
            % 那此时 LMP + Para.signI' .* Var.spread = Para.xpara' - Var.dual_mumax + Var.dual_mumin
           
            % Var.MX = (Para.xpara' - Var.dual_mumax + Var.dual_mumin) * Var.QX 
            if isfield(Setting,'new_version') && (Setting.new_version ~= 0)
                % 如果是new_version的话，整体要按照qabsmax去进行kIR的分摊
                Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Para.qmax(tempno) * Var.dual_mumax(tempno) ...
                                + Para.qmin(tempno) *  Var.dual_mumin(tempno) + Var.kIR * sum(Para.qabsmax(tempno))];
            else
                Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Para.qmax(tempno) * Var.dual_mumax(tempno) ...
                                + Para.qmin(tempno) *  Var.dual_mumin(tempno)];
            end
            % 如果是买家的话，买家出价 > 边际价格，此时达到qmax, 结算价格要略减，没问题
            Cons = [Cons, Var.RX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Var.MX(c,1)];      
        end 
    end

    
    % Calculate BB. 根据情况看是否要加上congestion revenue
    if Setting.topology == 0 || Setting.direct_solve == 1%如果不含有拓扑或者直接求解
        Cons = [Cons, sum(Var.MX) == Var.BB];
    else
        Var.congestion_revenue = sdpvar(1,1);
        Cons = [Cons, Var.congestion_revenue == sum(Var.dual_rhomax .* Para.BranchMax) + ...
            sum(Var.dual_rhomin .* Para.BranchMax)];
        Cons = [Cons, sum(Var.MX) - Var.congestion_revenue == Var.BB];
    end

    % Definition of RelaxBB, Relax
%      Cons = [Cons, Para.SW_max - Var.SW <= Var.RelaxSW];


    %% The relationship between RelaxIR and RelaxIC
    if isfield(Setting,'new_version') && (Setting.new_version ~= 0)
%         Cons = [Cons, Var.RelaxIR >= 0];
        for c = 1:Num.C
            tempno = Para.company_set(c).no;
            Cons = [Cons, sum(Para.qabsmax(tempno).* Para.xpara(tempno)) * Var.RelaxIR_C(c) + Var.RX(c) >= 0];
        end
        Cons = [Cons, Var.RelaxIR_C >=0];
        Cons = [Cons, Var.RelaxIR == sum(Var.RelaxIR_C)];
        
        Cons = [Cons, Var.RelaxIC >= 0];
%         Cons = [Cons, Var.RelaxIC >= Var.kIC/min(Para.xpara)];  %S^{IC} = 1/k
%         Cons = [Cons, Var.kIC >= 0];
        Cons = [Cons, Var.RelaxIC >= Var.kIC];
    end
    Cons = [Cons, Var.BB == -Var.RelaxBB];
end


%     Var.RelaxBB = sdpvar(1,1);
%     Var.RelaxIC = sdpvar(1,1);
%     Var.RelaxIR = sdpvar(1,1);
%     Var.RelaxSW = sdpvar(1,1);
%     Var.qip = sdpvar(Num.I,Num.P+1,'full');
% 
%     Var.mip = sdpvar(Num.I,Num.P+1,'full');
%     Var.rip = sdpvar(Num.I,Num.P+1,'full');
% %     Var.qip_mean = sdpvar(Num.I,Num.P,'full');
%     Var.qip_mean = sdpvar(Num.I, Num.P+2, 'full');
% 
% 
%     Cons = [];
%     
%     
%     % 这里是一些特殊的设计
%     
%     
%     % 直接对qipp的规则写死
%     % Cons = [Cons,Var.qip(1,:) == Para.Point(1,:)];
%     % Cons = [Cons,Var.qip(2,:) == Para.Point(1,:)-1];
%     
%     
%     if Setting.Qipp == 1 %应该必须要设计Qipp,否则的话qipp自由度太高
%         
%         %通过Qip来计算qip
%         %根据场景集来编号
%         Num.S = (Num.P+1) ^ Num.I;
%         Var.Qipp = sdpvar(Num.I,Num.S);
%         
%         
%         for i = 1:Num.I
%             Cons = [Cons, Var.Qipp(i,:) >= Para.qmin(i)];
%             Cons = [Cons, Var.Qipp(i,:) <= Para.qmax(i)];
%         end
%         
%         Cons = [Cons, sum(Var.Qipp,1) == 0];
%         Cons_num.Qbalance = length(Cons);
%         
%         % 现在可以再写一个拓扑约束. 拓扑约束必须要在Qipp=1的时候才成立
%         if isfield(Setting,'topology')
%             if Setting.topology == 1
%             % 写出一个GSDF. GSDF是什么, Num.Branch*Num.I, Num.I*Num.S
%                 Cons = [Cons, Para.GSDF * Var.Qipp <= repmat(Para.BranchMax,1,Num.S)];
%                 Cons = [Cons, Para.GSDF * Var.Qipp >= repmat(-Para.BranchMax,1,Num.S)];
%             end
%         end
%         
%         
%         scenario_prob = ones(1,Num.S);
%         for i = 1:Num.I
%             scenario_set(i,:) =  floor(mod(0:Num.S-1,(Num.P+1)^(Num.I- i +1))/((Num.P+1)^(Num.I-i)));
%             scenario_prob =  scenario_prob .* Para.Prob(i,scenario_set(i,:)+1);
%         end
%         
%         
%         
%         for i = 1:Num.I
%             for p = 1:Num.P+1
%                 Cons = [Cons, Var.qip(i,p) == sum(Var.Qipp(i,find(scenario_set(i,:)==p-1)) .* scenario_prob(find(scenario_set(i,:)==p-1)) / Para.Prob(i,p))];
%                 if abs(Para.Prob(i,p) - sum(scenario_prob(find(scenario_set(i,:)==p-1)))) > 1e-10
%                     fprintf('fail %d %d',i,p)
%                     disp(i)
%                 end 
%                 Cons_num.Qqlink(i,p) = length(Cons);
%             end 
%         end       
%         
% %         Para.Point(scenario_set * Num.I + repmat([1:Num.I]',1,Num.S));;
% %         % 如果将出清规则参数化,那么出清模型需要设置一个KKT条件.         
% %         if Setting.parameterize == 1
% %             Var.ramsey_price = sdpvar(Num.I,1); %设置的等效信息价格
% %             
% %             Cons = [Cons, Var.ramsey_price(Para.Gset) >=0];
% %             Cons = [Cons, Var.ramsey_price(Para.Dset) <=0];         
% %         end 
% 
%           % 这是本来双边交易时简单的写法
% %         Var.Qipp = sdpvar(Num.I,Num.P+1,Num.P+1);
% %         for p = 1:Num.P+1
% %             Cons = [Cons, Var.qip(1,p) == sum(reshape(Var.Qipp(1,p,:),1,Num.P+1) .* Para.Prob(2,:))];
% %             Cons = [Cons, Var.qip(2,p) == sum(reshape(Var.Qipp(2,:,p),1,Num.P+1) .* Para.Prob(1,:))];
% %         end 
% %         Cons = [Cons, Var.Qipp(1,:,:) + Var.Qipp(2,:,:) == 0];
% %         for i = 1:Num.I
% %             Cons = [Cons, Var.Qipp(i,:,:) >= Para.qmin(i)];
% %             Cons = [Cons, Var.Qipp(i,:,:) <= Para.qmax(i)];
% %         end
%     end 
%     if Setting.Qipp == 0        
%         if strcmp(Setting.appro.method,'inverse') % 如果采用伪逆的做法近似
%     
%     
%     % 直接把Qipp的把规则限制死
%     % for p1 = 1:Num.P+1 %买方，卖方
%     %     for p2= 1:Num.P+1
%     %         if p1>p2
%     %             Cons = [Cons, Var.Qipp(1,p1,p2) == 1];
%     %             Cons = [Cons, Var.Qipp(2,p1,p2) == -1];
%     %         elseif p1 == p2
%     %             Cons = [Cons, Var.Qipp(1,p1,p2) == 0.5];
%     %             Cons = [Cons, Var.Qipp(2,p1,p2) == -0.5];
%     %         else
%     %             Cons = [Cons, Var.Qipp(1,p1,p2) == 0];
%     %             Cons = [Cons, Var.Qipp(2,p1,p2) == 0];
%     %         end
%     %     end 
%     % end 
%     
%     
%     %社会福利. 只有这个社会福利约束最难达到. 这里把一开始的SW_max直接替换为0
%     Cons = [Cons, sum(sum(Para.Point .* Para.Prob .* Var.qip)) >= 0 - Var.RelaxSW - Setting.SWtolerance];
%     Cons_num.welafre = length(Cons);
%     % 这个地方感觉应该改一下
%     Cons = [Cons, (Var.qip(:,2:Num.P+1) + Var.qip(:,1:Num.P))/2 == Var.qip_mean(:,2:Num.P+1)]; %计算interval中间的qip_mean值.
%     Cons = [Cons, Var.qip(:,1) == Var.qip_mean(:,1)];
%     Cons = [Cons, Var.qip(:,Num.P+1) == Var.qip_mean(:,Num.P+2)];
%     Cons_num.qipmean = length(Cons);
%     Cons = [Cons, sum(sum(Para.Prob .* Var.mip)) >= -Var.RelaxBB];    %收支平衡约束
%     Cons_num.mbalance = length(Cons);
%     Cons = [Cons, sum(sum(Para.Prob .* Var.qip)) == 0];               %物品平衡约束
%     Cons_num.qbalance = length(Cons);
% 
%     for i = 1:Num.I
%         % 容量上下限约束
%         Cons = [Cons, Var.qip(i,:) >= Para.qmin(i) ];
%         Cons = [Cons, Var.qip(i,:) <= Para.qmax(i) ];
%         for p = 1:Num.P+1
%             if p > 1
%                 Cons = [Cons, Var.qip(i,p) >= Var.qip(i,p-1)]; %单调递增约束
%             end 
%             if Para.type(i) == 1 %如果是买方  
% %                  Cons = [Cons, Var.rip(i,p) == sum(sum((Var.qip_mean(i,1:p-1) - Var.RelaxIC) .* Para.Interval_len(i,1:p-1))) - Var.RelaxIR];
%                  % 有的时候IR不一定是参数最差的时候取值
%                  Cons = [Cons, Var.rip(i,p) == sum(sum((Var.qip_mean(i,1:p) - Var.RelaxIC) .* Para.Interval_len(i,1:p))) - Var.RelaxIR];
%                  Cons = [Cons, Var.mip(i,p) == Para.Point(i,p) * Var.qip(i,p) - Var.rip(i,p)];
%             end
%             if Para.type(i) == 2 %如果是卖方
% %                  Cons = [Cons, Var.rip(i,p) == sum(sum((-Var.qip_mean(i,p:Num.P) - Var.RelaxIC) .* Para.Interval_len(i,p:Num.P))) - Var.RelaxIR];
%                  Cons = [Cons, Var.rip(i,p) == sum(sum((-Var.qip_mean(i,p+1:Num.P+2) - Var.RelaxIC) .* Para.Interval_len(i,p+1:Num.P+2))) - Var.RelaxIR];
%                  Cons = [Cons, Var.mip(i,p) == Para.Point(i,p) * Var.qip(i,p) - Var.rip(i,p)];
%             end
%         end
%     end
%     if isfield(Setting,'topology')
%         if Setting.topology == 1
%             % 写出一个GSDF. GSDF是什么, Num.Branch*Num.I, Num.I*Num.S
%             Cons = [Cons, Para.GSDF * Var.qip <= repmat(Para.BranchMax,1,Num.P+1)];
%             Cons = [Cons, Para.GSDF * Var.qip >= repmat(-Para.BranchMax,1,Num.P+1)];
%         end
%    end


    
