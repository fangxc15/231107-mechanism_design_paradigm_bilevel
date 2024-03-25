function  Result = optimization_primal(Num,Para,Setting, varargin)
    % 原问题的最优化
    % 用LMP结算时，到底要不要扣除spread是有讲究的
    if isfield(Setting,'falsebid_volume') %增加falsebid
        bidPara = Para;
        bidPara.xpara = bidPara.xpara + Setting.falsebid_volume(:)';
        bidPara = process_Para(bidPara,Num,Setting);
    else
        bidPara = Para;
    end
    
        
    if length(varargin) >= 1
        input_spread = cell2mat(varargin(1));
    end
    if length(varargin) >= 2
        kIR = cell2mat(varargin(2));
    end
    if length(varargin) >= 3
        kIC = cell2mat(varargin(3));
    end
    % Variables for lower level
    for c = 1:Num.C 
        Var.scenario(c).QX = sdpvar(Num.I,Num.P,'full');
        Var.scenario(c).QXmean = sdpvar(Num.I,Num.P,'full');
    end 
    Var.QX = sdpvar(Num.I,1);
%     Var.MX = sdpvar(Num.C,1); 
%     Var.RX = sdpvar(Num.C,1); 

    Var.RelaxBB = sdpvar(1,1);
    Var.RelaxSW = sdpvar(1,1);
    Var.RelaxIC = sdpvar(1,1);
    Var.RelaxIR = sdpvar(1,1);
%     Var.BB = sdpvar(1,1);
    Var.obj = sdpvar(1,1);
    Var.SW = sdpvar(1,1);
    Var.Gcost = sdpvar(1,1);
    Var.Dutil = sdpvar(1,1);
    Var.bidSW = sdpvar(1,1);
    Var.bidGcost = sdpvar(1,1);
    Var.bidDutil = sdpvar(1,1);
%     Var.Grevenue = sdpvar(1,1);
%     Var.Drevenue = sdpvar(1,1);
    
    Cons = [];
    
    % Base scenario
    Cons = [Cons, Var.QX - bidPara.qmin' >=0]; Cons_no.dual_mumin = length(Cons);
    Cons = [Cons, bidPara.qmax' - Var.QX >=0]; Cons_no.dual_mumax = length(Cons);
    Cons = [Cons, sum(Var.QX) <= 0];        Cons_no.dualprice = length(Cons); %这个其实表示了需求一定小于供给
    
    if Setting.topology == 1
        Cons = [Cons, bidPara.BranchMax + bidPara.GSDF * Var.QX >=0];   Cons_no.dual_rhomax = length(Cons);
        Cons = [Cons, - bidPara.GSDF * Var.QX + bidPara.BranchMax >=0]; Cons_no.dual_rhomin = length(Cons);
    end

    
    % Define welfare
    if Setting.riskSW == 0
        Cons = [Cons, Para.xpara(Para.Gset) * Var.QX(Para.Gset) == Var.Gcost];
        Cons = [Cons, Para.xpara(Para.Dset) * Var.QX(Para.Dset) == Var.Dutil]; 
        Cons = [Cons, Var.Gcost + Var.Dutil == Var.SW];
        
        Cons = [Cons, bidPara.xpara(bidPara.Gset) * Var.QX(bidPara.Gset) == Var.bidGcost];
        Cons = [Cons, bidPara.xpara(bidPara.Dset) * Var.QX(bidPara.Dset) == Var.bidDutil];    
        Cons = [Cons, Var.bidGcost + Var.bidDutil == Var.bidSW];

        Cons = [Cons, (bidPara.xpara - bidPara.sign.*input_spread) * Var.QX == Var.obj]; % 在计算bidding的时候，这里是减号。卖家-1，相当于卖家向上加
    else
        Var.tempSW = sdpvar(Num.C, Num.P,'full');
        Var.tempGcost = sdpvar(Num.C, Num.P,'full');
        Var.tempDutil = sdpvar(Num.C, Num.P,'full');
        Var.baseSW = sdpvar(1,1);
        Var.baseGcost = sdpvar(1,1);
        Var.baseDutil = sdpvar(1,1);
        
        Var.tempbidSW = sdpvar(Num.C, Num.P,'full');
        Var.tempbidGcost = sdpvar(Num.C, Num.P,'full');
        Var.tempbidDutil = sdpvar(Num.C, Num.P,'full');
        Var.basebidSW = sdpvar(1,1);
        Var.basebidGcost = sdpvar(1,1);
        Var.basebidDutil = sdpvar(1,1);
        
        Var.baseobj = sdpvar(1,1);
        Var.tempobj = sdpvar(Num.C, Num.P, 'full');
        
        for c = 1:Num.C
            if bidPara.Cstrategic(c) == 0
                Cons = [Cons, Var.scenario(c).QX == 0];
            end
            for p = 1:Num.P
                tempno = bidPara.company_set(c).no;
                temppara = bidPara.xpara;
                temppara(tempno) = bidPara.Point(tempno, p+1);
                
                temp_realpara = Para.xpara;
                temp_realpara(tempno) = Para.Point(tempno,p+1);
                
                Cons = [Cons, temp_realpara * Var.scenario(c).QX(:,p) == Var.tempSW(c,p)];
                Cons = [Cons, temp_realpara(Para.Gset) * Var.scenario(c).QX(Para.Gset,p) == Var.tempGcost(c,p)];
                Cons = [Cons, temp_realpara(Para.Dset) * Var.scenario(c).QX(Para.Dset,p) == Var.tempDutil(c,p)];
                
                Cons = [Cons, temppara * Var.scenario(c).QX(:,p) == Var.tempbidSW(c,p)];
                Cons = [Cons, temppara(bidPara.Gset) * Var.scenario(c).QX(bidPara.Gset,p) == Var.tempbidGcost(c,p)];
                Cons = [Cons, temppara(bidPara.Dset) * Var.scenario(c).QX(bidPara.Dset,p) == Var.tempbidDutil(c,p)];

                
                Cons = [Cons, Var.tempobj(c,p) == (temppara  - bidPara.signI.* input_spread) * Var.scenario(c).QX(:,p)];
            end
        end
        Cons = [Cons, Para.xpara(Para.Gset) * Var.QX(Para.Gset) == Var.baseGcost];
        Cons = [Cons, Para.xpara(Para.Dset) * Var.QX(Para.Dset) == Var.baseDutil];
        Cons = [Cons, Var.baseSW == Para.xpara * Var.QX];
        
        Cons = [Cons, bidPara.xpara(bidPara.Gset) * Var.QX(bidPara.Gset) == Var.basebidGcost];
        Cons = [Cons, bidPara.xpara(bidPara.Dset) * Var.QX(bidPara.Dset) == Var.basebidDutil];
        Cons = [Cons, Var.basebidSW == bidPara.xpara * Var.QX];
        
        Cons = [Cons, Var.baseobj == (bidPara.xpara  - bidPara.signI.* input_spread) * Var.QX];
        
        Cons = [Cons, Var.SW == Para.prob(1,1) * Var.baseSW + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempSW))];
        Cons = [Cons, Var.Gcost == Para.prob(1,1) * Var.baseGcost + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempGcost))];
        Cons = [Cons, Var.Dutil == Para.prob(1,1) * Var.baseDutil + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempDutil))];
        
        Cons = [Cons, Var.bidSW == bidPara.prob(1,1) * Var.basebidSW + sum(sum(bidPara.prob(:,2:Num.P+1) .* Var.tempbidSW))];
        Cons = [Cons, Var.bidGcost == bidPara.prob(1,1) * Var.basebidGcost + sum(sum(bidPara.prob(:,2:Num.P+1) .* Var.tempbidGcost))];
        Cons = [Cons, Var.bidDutil == bidPara.prob(1,1) * Var.basebidDutil + sum(sum(bidPara.prob(:,2:Num.P+1) .* Var.tempbidDutil))];

        
        Cons = [Cons, Var.obj == Var.baseobj + sum(sum(Var.tempobj))];
    end
    % Calculate QX mean
    for c = 1:Num.C 
        Cons = [Cons, Var.scenario(c).QXmean(:,1) == (Var.scenario(c).QX(:,1) + Var.QX)/2];
        for p = 2:Num.P
            Cons = [Cons, Var.scenario(c).QXmean(:,p) == (Var.scenario(c).QX(:,p) + Var.scenario(c).QX(:,p-1))/2];
        end
    end 
    
    
    % Extended scenario
    for c = 1:Num.C
        if bidPara.Cstrategic(c) == 1
            Cons = [Cons, sum(Var.scenario(c).QX,1) <= 0];                      Cons_no.scenario(c).dualprice = length(Cons);
            Cons = [Cons, Var.scenario(c).QX - repmat(bidPara.qmin',1,Num.P) >=0]; Cons_no.scenario(c).dual_mumin = length(Cons);
            Cons = [Cons, repmat(bidPara.qmax',1,Num.P) - Var.scenario(c).QX >=0]; Cons_no.scenario(c).dual_mumax = length(Cons);
            for p = 1:Num.P
                temppara = bidPara.xpara;
                tempno = bidPara.company_set(c).no;
                temppara(tempno) = bidPara.Point(tempno, p+1);
                if Setting.topology == 1
                    Cons = [Cons, bidPara.BranchMax + bidPara.GSDF * Var.scenario(c).QX(:,p) >=0]; Cons_no.scenario(c).dual_rhomax(p) = length(Cons);
                    Cons = [Cons, bidPara.BranchMax - bidPara.GSDF * Var.scenario(c).QX(:,p) >=0]; Cons_no.scenario(c).dual_rhomin(p) = length(Cons);
                end
            end
        else
            Cons = [Cons, Var.scenario(c).QX == 0];
        end
    end
    
%     % Calculate MX
%     for c = 1:Num.C
%         tempno = bidPara.company_set(c).no;
%                
%         if bidPara.Cstrategic(c) == 1
%             if length(tempno) > 1 % a company with multiple units
%                 temp_len = sum(sum(bidPara.Interval_len(tempno,:) .^2,1).^0.5);
%                 Cons = [Cons, Var.RX(c,1) == sum(sum(bidPara.Interval_len(tempno,:) .* ...
%                     Var.scenario(c).QXmean(tempno,:))) - Var.RelaxIC * temp_len - Var.RelaxIR];
%                 Cons = [Cons, Var.MX(c,1) == bidPara.xpara(tempno) * Var.QX(tempno) - Var.RX(c,1)];
%     %             Cons = [Cons, Var.MX(c,1) == bidPara.xpara(tempno) * Var.QX(tempno) - sum(sum(bidPara.Interval_len(tempno,:) .* ...
%     %                 Var.scenario(c).QXmean(tempno,:))) + bidPara.sign(c) * Var.RelaxIC * temp_len + Var.RelaxIR];
%             else % a single unit
%                 Cons = [Cons, Var.RX(c,1) == sum(bidPara.Interval_len(tempno,:) .* ...
%                     (Var.scenario(c).QXmean(tempno,:) - bidPara.sign(c) * Var.RelaxIC)) - Var.RelaxIR];
%                 Cons = [Cons, Var.MX(c,1) == bidPara.xpara(tempno) * Var.QX(tempno) - Var.RX(c,1)];
%     %             Cons = [Cons, Var.MX(c,1) == bidPara.xpara(tempno) * Var.QX(tempno) - sum(bidPara.Interval_len(tempno,:) .* ...
%     %                 (Var.scenario(c).QXmean(tempno,:) - bidPara.sign(c) * Var.RelaxIC)) + Var.RelaxIR];
%             end
%             
%         else
%             Cons = [Cons, Var.MX(c,1) == bidPara.xpara(tempno) * Var.QX(tempno) - bidPara.qmax(tempno) * Var.dual_mumax(tempno) ...
%                                 + bidPara.qmin(tempno) *  Var.dual_mumin(tempno)];
%             % 如果是买家的话，买家出价 > 边际价格，此时达到qmax, 结算价格要略减，没问题
%             Cons = [Cons, Var.RX(c,1) == bidPara.xpara(tempno) * Var.QX(tempno) - Var.MX(c,1)];      
%         end 
%     end

    
%     % Calculate BB
%     if Setting.topology == 0 || Setting.direct_solve == 1%如果不含有拓扑或者直接求解
%         Cons = [Cons, sum(Var.MX) == Var.BB];
%     else
%         Var.congestion_revenue = sdpvar(1,1);
%         Cons = [Cons, Var.congestion_revenue == sum(Var.dual_rhomax .* bidPara.BranchMax) + ...
%             sum(Var.dual_rhomin .* bidPara.BranchMax)];
%         Cons = [Cons, sum(Var.MX) - Var.congestion_revenue == Var.BB];
%     end

%     Cons = [Cons, Var.BB == -Var.RelaxBB];
    obj = -Var.obj;
    
    ops = sdpsettings('solver','gurobi');
    ops.gurobi.MIPFocus = 0;      %默认值0,试图在最优值和可行解之间取得平衡。1：以可行解为目标；2：以最优解为目标；3：最优边界为目标。
    ops.gurobi.MIPGap = 0.001;   %默认值0.0001,从某Gap就停止，相对值
    ops.gurobi.MIRCuts = 2;       %Aggressive (2), Conservative (1), Automatic (-1), or None (0)
    ops.gurobi.Cuts = 3;          %0 to shut off cuts, 1 for moderate cut generation, 2 for aggressive cut generation, and 3 for very aggressive cut generation
    ops.gurobi.ImproveStartTime = 300;
    ops.gurobi.ImproveStartGap = 0.01;
    ops.gurobi.TimeLimit = 600;
    ops.cplex.MIPGap = 0.001; 
    solution = optimize(Cons,obj,ops);
    
    Result.obj = value(obj);
    
    % Get spread
    Result.spread = input_spread;
    Result.unified_spread = abs(input_spread(1));

    % 原变量和对偶变量
    Result.QX = value(Var.QX);
    for c = 1:Num.C 
        Result.scenario(c).QX = value(Var.scenario(c).QX);
        Result.scenario(c).QXmean = value(Var.scenario(c).QXmean);
    end 
    if Setting.riskSW == 1
        Result.qtotal =  sum(abs(Result.QX)) * Para.prob(1,1);
        for c = 1:Num.C
            for p = 1:Num.P
                Result.qtotal = Result.qtotal + Para.prob(c,p+1) * sum(abs(Result.scenario(c).QX(:,p)));
            end
        end  
    else
        Result.qtotal = sum(abs(Result.QX));
    end
    
    Result.dualprice = dual(Cons(Cons_no.dualprice));
    Result.dual_mumax = dual(Cons(Cons_no.dual_mumax));
    Result.dual_mumin = dual(Cons(Cons_no.dual_mumin));
    Result.dual_mumax_bin = Result.dual_mumax > 0;
    Result.dual_mumin_bin = Result.dual_mumin > 0;
    
    for c = 1:Num.C
        if bidPara.Cstrategic(c) == 1
            Result.scenario(c).dualprice = reshape(dual(Cons(Cons_no.scenario(c).dualprice)),1,Num.P);
            Result.scenario(c).dual_mumax = reshape(dual(Cons(Cons_no.scenario(c).dual_mumax)),Num.I,Num.P);
            Result.scenario(c).dual_mumin = reshape(dual(Cons(Cons_no.scenario(c).dual_mumin)),Num.I,Num.P);
            Result.scenario(c).dual_mumax_bin = Result.scenario(c).dual_mumax >0;
            Result.scenario(c).dual_mumin_bin = Result.scenario(c).dual_mumin >0;
        end
    end

    if Setting.topology == 1
        Result.dual_rhomax = dual(Cons(Cons_no.dual_rhomax));
        Result.dual_rhomin = dual(Cons(Cons_no.dual_rhomin));
        Result.dual_rhomax_bin = Result.dual_rhomax > 0;
        Result.dual_rhomin_bin = Result.dual_rhomin > 0;
        Result.Pl  =  -bidPara.GSDF * Result.QX;
        Result.LMP = Result.dualprice + bidPara.GSDF' * (Result.dual_rhomin - Result.dual_rhomax);
        for c = 1:Num.C
            if bidPara.Cstrategic(c) == 1  
                for p = 1:Num.P
                    Result.scenario(c).dual_rhomax(:,p) = dual(Cons(Cons_no.scenario(c).dual_rhomax(p)));
                    Result.scenario(c).dual_rhomin(:,p) = dual(Cons(Cons_no.scenario(c).dual_rhomin(p)));
                    Result.scenario(c).dual_rhomax_bin(:,p) = Result.scenario(c).dual_rhomax(:,p) > 0;
                    Result.scenario(c).dual_rhomin_bin(:,p) = Result.scenario(c).dual_rhomin(:,p) > 0;
%                     bidPara.GSDF * Var.scenario(c).QX(:,p) == -Var.scenario(c).Pl(:,p)              
                    Result.scenario(c).Pl(:,p)  =  -bidPara.GSDF * Result.scenario(c).QX(:,p);
%                     Result.scenario(c).LMP =  Result.dualprice + bidPara.GSDF' * (Result.dual_rhomin - Result.dual_rhomax);

                    Result.scenario(c).LMP(:,p) = - bidPara.GSDF' * (Result.scenario(c).dual_rhomax(:,p) - ...
                                        Result.scenario(c).dual_rhomin(:,p)) + repmat(Result.scenario(c).dualprice(1,p),Num.I,1);
                end
            end
        end
    else
        Result.LMP = repmat(Result.dualprice,Num.I,1);
        for c = 1:Num.C
            if bidPara.Cstrategic(c) == 1  
                for p = 1:Num.P
                    Result.scenario(c).LMP(:,p) = repmat(Result.scenario(c).dualprice(1,p),Num.I,1);
                end
            end
        end
    end
  
%         Cons = [Cons, Var.SW == bidPara.prob(1,1) * Var.baseSW + sum(sum(bidPara.prob(:,2:Num.P+1) .* Var.tempSW))];
%         Cons = [Cons, Var.Gcost == bidPara.prob(1,1) * Var.baseGcost + sum(sum(bidPara.prob(:,2:Num.P+1) .* Var.tempGcost))];
%         Cons = [Cons, Var.Dutil == bidPara.prob(1,1) * Var.baseDutil + sum(sum(bidPara.prob(:,2:Num.P+1) .* Var.tempDutil))];
%         Cons = [Cons, Var.Gcost + Var.Dutil == Var.sumSW];

    Result.clearing_obj = value(Var.obj);
    if Setting.riskSW == 0
        Result.SW = value(Var.SW);
        Result.Gcost = value(Var.Gcost);
        Result.Dutil = value(Var.Dutil);
        Result.bidSW = value(Var.bidSW);
        Result.bidGcost = value(Var.bidGcost);
        Result.bidDutil = value(Var.bidDutil);
    elseif Setting.riskSW == 1
        Result.baseobj = value(Var.baseobj);
        Result.baseSW = value(Var.baseSW);
        Result.baseGcost = value(Var.baseGcost);
        Result.baseDutil = value(Var.baseDutil);
        Result.tempobj = value(Var.tempobj);
        Result.tempSW = value(Var.tempSW);
        Result.tempGcost = value(Var.tempGcost);
        Result.tempDutil = value(Var.tempDutil);
        Result.SW = value(Var.SW);
        Result.Gcost = value(Var.Gcost);
        Result.Dutil = value(Var.Dutil);
        
        Result.basebidSW = value(Var.basebidSW);
        Result.basebidGcost = value(Var.basebidGcost);
        Result.basebidDutil = value(Var.basebidDutil);
        Result.tempbidSW = value(Var.tempbidSW);
        Result.tempbidGcost = value(Var.tempbidGcost);
        Result.tempbidDutil = value(Var.tempbidDutil);
        Result.bidSW = value(Var.bidSW);
        Result.bidGcost = value(Var.bidGcost);
        Result.bidDutil = value(Var.bidDutil);
%         Result.SW = bidPara.prob(1,1) * Result.baseSW + sum(sum(bidPara.prob(:,2:Num.P+1) .* Result.tempSW));
%         Result.Gcost = bidPara.prob(1,1) * Result.baseGcost + sum(sum(bidPara.prob(:,2:Num.P+1) .* Result.tempGcost));
%         Result.Dutil = bidPara.prob(1,1) * Result.baseDutil + sum(sum(bidPara.prob(:,2:Num.P+1) .* Result.tempDutil));      
    end
    Result.welfare = Result.SW;
    
    
    if isfield(Setting,'new_version') && (Setting.new_version ~=0)
        % 在new_version里，需要考虑传入的kIC和kIR项
        if exist('kIC')
            Result.kIC = kIC;
        else
            Result.kIC = 0;
        end
        if exist('kIR')
            Result.kIR = kIR;
        else
            Result.kIR = 0;
        end
        
        Result = calculate_MX(Result, Setting, Num, Para);
        Result = calculate_Relax(Result, Setting, Num, Para);
%         if Result.unified_spread == 0
%             bidPara.SW_max = Result.SW;
%             Result.RelaxSW = 0;
%         else
%             Result.RelaxSW = bidPara.SW_max - Result.SW;
%         end
        
%         Result.RelaxIC = Result.kIC/min(bidPara.xpara);
%         Result.RelaxIC = max(0,Result.RelaxIC);
        
%         for c = 1:Num.C
%             tempno = bidPara.company_set(c).no;
%             Result.RelaxIR_C(c) = max(0,-Result.RX(c)/ sum(bidPara.qabsmax(tempno).* bidPara.xpara(tempno)));
%         end
%         Result.RelaxIR = sum(Result.RelaxIR_C);
            
% %             if isfield(Result,'RelaxIR')
% %                 Result.RelaxIR = max(Result.RelaxIR, -Result.RX(c)/ sum(bidPara.qabsmax(tempno).* bidPara.xpara(tempno)));
% %             else
% %                 Result.RelaxIR = -Result.RX(c)/ sum(bidPara.qabsmax(tempno).* bidPara.xpara(tempno));
% %             end
               
    else
        Result.RelaxIC = 0;
        Result.RelaxIR = 0;
        Result.RelaxSW = -Result.SW;
        Result = calculate_MX(Result, Setting, Num, Para);
    end



    
    
%     if Setting.direct_solve == 0
%         if isfield(Setting,'settle_spread') && Setting.settle_spread == 0
%             Result.LMPI = Result.LMP; %否则按照spread结算
%         else
%             % sign的话卖家是-1
%             Result.LMPI = (Result.LMP + bidPara.signI' .* Result.spread(:));  % 如果没有settle spread，其实默认按照LMP减掉申报的spread
%         end
%         Result.LMPsettle = Result.LMPI .* Result.QX;
%         Result.LMPsettleC = zeros(Num.C,1);
%         for c = 1:Num.C
%             Result.LMPsettleC(c) = sum(Result.LMPsettle(bidPara.company_set(c).no));
%         end
%         Result.Settle_diff = Result.MX - Result.LMPsettleC;
%     end
% 
%     if Setting.topology == 1 && Setting.direct_solve == 0
%         Result.LMPbalance = sum(Result.LMPsettleC) - Result.congestion_revenue;
%     end
    
%     Result.QXPara = zeros(Num.C,1);
%     for c = 1:Num.C %这个可以补充进来, 但这里指的是基础状态下的生产成本
%         tempno = bidPara.company_set(c).no;
%         Result.QXPara(c) = bidPara.xpara(tempno) * Result.QX(tempno);       
%     end
%     Result.producer_welfare = sum(Result.RX(bidPara.Gset_company));
%     Result.consumer_welfare = sum(Result.RX(bidPara.Dset_company));
%     if Setting.direct_solve == 0
%         Result.LMPsettleC_RX = Result.QXPara - Result.LMPsettleC;
%     end
end
%     if Setting.riskSW == 1
%         Result.qtotal =  sum(abs(Result.QX)) * bidPara.prob(1,1);
%         for c = 1:Num.C
%             for p = 1:Num.P
%                 Result.qtotal = Result.qtotal + bidPara.prob(c,p+1) * sum(abs(Result.scenario(c).QX(:,p)));
%             end
%         end  
%     else
%         Result.qtotal = sum(abs(Result.QX));
%     end
%     Result.surplus = Result.BB;


%     if Setting.direct_solve == 0
%         Result.LMPI = (Result.LMP + bidPara.signI' .* Result.spread);
%         Result.LMPsettle = Result.LMPI .* Result.QX;
%         Result.LMPsettleC = zeros(Num.C,1);
%         for c = 1:Num.C
%             Result.LMPsettleC(c) = sum(Result.LMPsettle(bidPara.company_set(c).no));
%         end
%         Result.Settle_diff = Result.MX - Result.LMPsettleC;
%     end
    
%     if Setting.topology == 1
%         Result.Pl = bidPara.GSDF * Result.QX;
%     end
    
