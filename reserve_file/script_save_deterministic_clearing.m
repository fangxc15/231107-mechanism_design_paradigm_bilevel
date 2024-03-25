% 在这里，我甚至无需去设计Q(X)单调递增的约束, 因为它肯定是单调递增的.
% 要注意互补松弛条件写的时候Para.bigM的数值问题


Var.spread = sdpvar(Num.I,1); %这个spread可以

Var.MX = sdpvar(Num.C,1); 
Var.RelaxBB = sdpvar(1,1);
Var.RelaxSW = sdpvar(1,1);
Var.RelaxIC = sdpvar(1,1);
Var.RelaxIR = sdpvar(1,1);
Var.BB = sdpvar(1,1);
Var.SW = sdpvar(1,1);

% Variables for lower level
for c = 1:Num.C 
    Var.scenario(c).QX = sdpvar(Num.I,Num.P,'full');
    Var.scenario(c).QXmean = sdpvar(Num.I,Num.P,'full');
end 
Var.QX = sdpvar(Num.I,1);

Cons = [];

% Constraints on spread design
Cons = [Cons, Var.spread >= 0];
if Setting.unified_spread == 1    
    Var.unified_spread = sdpvar(1,1);
    Cons = [Cons, Var.spread == Var.unified_spread];
end 


% Calculate Welfare
if Setting.riskSW == 0
    Cons = [Cons, Para.xpara * Var.QX == Var.SW];
else
    Var.tempSW = sdpvar(Num.C, Num.P,'full');
    Var.baseSW = sdpvar(1,1);
    for c = 1:Num.C
        for p = 1:Num.P
            temppara = Para.xpara;
            tempno = Para.company_set(c).no;
            temppara(tempno) = Para.Point(tempno, p+1);
            Cons = [Cons, temppara * Var.scenario(c).QX(:,p) == Var.tempSW(c,p)];
        end
    end
    Cons = [Cons, Var.baseSW == Para.xpara * Var.QX];
    Cons = [Cons, Var.SW == Para.prob(1,1) * Var.baseSW + sum(sum(Para.prob(:,2:Num.P+1) .* Var.tempSW))];
end

% Calculate QX mean
for c = 1:Num.C 
    Cons = [Cons, Var.scenario(c).QXmean(:,1) == (Var.scenario(c).QX(:,1) + Var.QX)/2];
    for p = 2:Num.P
        Cons = [Cons, Var.scenario(c).QXmean(:,p) == (Var.scenario(c).QX(:,p) + Var.scenario(c).QX(:,p-1))/2];
    end
end 

% Calculate MX
for c = 1:Num.C
    tempno = Para.company_set(c).no;
    if length(tempno) > 1 % a company with multiple units
        temp_len = sum(sum(Para.Interval_len(tempno,:) .^2).^0.5);
        Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - sum(sum(Para.Interval_len(tempno,:) .* ...
            Var.scenario(c).QXmean(tempno,:))) + Para.sign(c) * Var.RelaxIC * temp_len + Var.RelaxIR];
    else % a single unit
        Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - sum(Para.Interval_len(tempno,:) .* ...
            (Var.scenario(c).QXmean(tempno,:) - Para.sign(c) * Var.RelaxIC)) + Var.RelaxIR];
    end
end

% Calculate BB
Cons = [Cons, sum(Var.MX) == Var.BB];

% Definition of RelaxBB, Relax
Cons = [Cons, Para.SW_max - Var.SW == Var.RelaxSW];
Cons = [Cons, Var.BB == -Var.RelaxBB];

% Tempobj
obj = Var.RelaxSW;
Cons = [Cons, Var.RelaxIC <=0, Var.RelaxIR <= 0,Var.RelaxBB <=0];


%% Direct Design of Lower Level Model
if Setting.direct_solve == 1
    % Lower Level Model Design
    obj_base = (Para.xpara  - Para.sign .* Var.spread')* Var.QX; % Clearing for one case
    obj_inner = obj_base;
    Cons_inner = [];
    Cons_inner = [Cons_inner, sum(Var.QX) == 0];
    Cons_inner = [Cons_inner, Var.QX <= Para.qmax', Var.QX >= Para.qmin'];
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

    
    
    

%% KKT conditions for lower level model
elseif Setting.direct_solve == 0
    Var.dualprice = sdpvar(1,1);
    Var.dual_mumax = sdpvar(Num.I,1);
    Var.dual_mumin = sdpvar(Num.I,1);
    Var.dual_mumax_bin = binvar(Num.I,1);
    Var.dual_mumin_bin = binvar(Num.I,1);
    if Setting.topology == 1
        Var.dual_rhomax = sdpvar(Num.Branch,1);
        Var.dual_rhomin = sdpvar(Num.Branch,1);
        Var.dual_rhomax_bin = binvar(Num.Branch,1);
        Var.dual_rhomin_bin = binvar(Num.Branch,1);
    end 
    Cons_inner = []
    Cons_inner = [Cons_inner, sum(Var.QX) == 0];
    % for i = 1:Num.I
    Cons_inner = F_Complementarity_Yal(Cons_inner, Var.QX - Para.qmin', Var.dual_mumin, Var.dual_mumin_bin);
    Cons_inner = F_Complementarity_Yal(Cons_inner, Para.qmax' - Var.QX, Var.dual_mumax, Var.dual_mumax_bin);
    %目前暂时的不含有拓扑的出清,Lagrange条件
    Cons_inner = [Cons_inner, Para.sign' .* Var.spread - Para.xpara' + Var.dual_mumax - Var.dual_mumin + Var.dualprice == 0];


    for c = 1:Num.C

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
        end
        Cons_inner = [Cons_inner, sum(Var.scenario(c).QX,1) == 0];
        Cons_inner = F_Complementarity_Yal(Cons_inner, Var.scenario(c).QX - repmat(Para.qmin',1,Num.P), Var.scenario(c).dual_mumin, Var.scenario(c).dual_mumin_bin);
        Cons_inner = F_Complementarity_Yal(Cons_inner, repmat(Para.qmax',1,Num.P) - Var.scenario(c).QX, Var.scenario(c).dual_mumax, Var.scenario(c).dual_mumax_bin);
        for p = 1:Num.P
            temppara = Para.xpara;
            tempno = Para.company_set(c).no;
            temppara(tempno) = Para.Point(tempno, p+1);
            Cons_inner = [Cons_inner, Para.sign' .* Var.spread - temppara' + Var.scenario(c).dual_mumax(:,p) ... 
                                            - Var.scenario(c).dual_mumin(:,p) + repmat(Var.scenario(c).dualprice(1,p),Num.I,1) == 0];
        end
    end 
    Cons = [Cons, Cons_inner];
    ops = sdpsettings('solver','gurobi');

    solution = optimize(Cons,obj,ops);
    

    
    for c = 1:Num.C
        Result.scenario(c).dual_mumin = value(Var.scenario(c).dual_mumin);
        Result.scenario(c).dual_mumax = value(Var.scenario(c).dual_mumax);
        Result.scenario(c).dualprice = value(Var.scenario(c).dualprice);
    end
    Result.dual_mumin = value(Var.dual_mumin);
    Result.dual_mumax = value(Var.dual_mumax);
    Result.dualprice = value(Var.dualprice);
end

Result.obj = value(obj);
Result.QX = value(Var.QX);
Result.MX = value(Var.MX);
Result.BB = value(Var.BB);
Result.SW = value(Var.SW);
if Setting.riskSW == 1
    Result.tempSW = value(Var.tempSW);
    Result.baseSW = value(Var.baseSW);
end
Result.RelaxBB = value(Var.RelaxBB);
Result.RelaxSW = value(Var.RelaxSW);
Result.RelaxIR = value(Var.RelaxIR);
Result.RelaxIC = value(Var.RelaxIC);
if Setting.unified_spread == 1    
    Result.unified_spread = value(Var.unified_spread);
end
Result.spread = value(Var.spread);

for c = 1:Num.C 
    Result.scenario(c).QX = value(Var.scenario(c).QX);
    Result.scenario(c).QXmean = value(Var.scenario(c).QXmean);
end 
    





%%
% c = 2;
% p = 1;
% temppara = Para.xpara;
% tempno = Para.company_set(c).no;
% temppara(tempno) = Para.Point(tempno, p+1);
% value(Para.sign' .* Var.spread - temppara' + Var.scenario(c).dual_mumax(:,p) ... 
%                                             - Var.scenario(c).dual_mumin(:,p) + Var.scenario(c).dualprice(1,p))
%                                         
% [value(Para.sign' .* Var.spread) -temppara' value(Var.scenario(c).dual_mumax(:,p)) ...
%     -value(Var.scenario(c).dual_mumin(:,p)) repmat(value(Var.scenario(c).dualprice(1,p)),Num.C,1)]
% 
% 
% 
% % 0.12   -0.60    0    0   0.48
% % -0.12  -0.44  0.08   0   0.48
% 
% % 0.12   -0.60    0   -0.08   0.56
% % -0.12  -0.44    0     0     0.56
% 
% 
% value(repmat(Para.qmax',1,Num.P) - Var.scenario(c).QX)
% value(Var.scenario(c).dual_mumax)
% value(Var.scenario(c).dual_mumax_bin)
% 
% value(Var.scenario(c).QX - repmat(Para.qmin',1,Num.P))
% value(Var.scenario(c).dual_mumin)
% value(Var.scenario(c).dual_mumin_bin)
% 
% value(Var.scenario(c).dual_mumin) == 0
% 1-value(Var.scenario(c).dual_mumin_bin)

