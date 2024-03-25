function Result = calculate_MX(Result, Setting, Num, Para)
    %% Calculate MX
    Result.bidRX = zeros(Num.C,1);
    Result.RX = zeros(Num.C,1);
    Result.MX = zeros(Num.C,1);
    
    if isfield(Setting,'falsebid_volume') %增加falsebid
        bidPara = Para;
        bidPara.xpara = bidPara.xpara + Setting.falsebid_volume;
        bidPara = process_Para(bidPara,Num,Setting);
    else
        bidPara = Para;
    end
    
    for c = 1:Num.C
        tempno = bidPara.company_set(c).no;
        
%         if Para.Cstrategic(c) == 1
%             if length(tempno) > 1 % a company with multiple units
%                 temp_len = sum(sum(Para.Interval_len(tempno,:) .^2,1).^0.5);
%                 Cons = [Cons, Var.RX(c,1) == sum(sum(Para.Interval_len(tempno,:) .* ...
%                     Var.scenario(c).QXmean(tempno,:))) - Var.RelaxIC * temp_len - Var.RelaxIR];
%                 Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Var.RX(c,1)];
%             else % a single unit
%                 Cons = [Cons, Var.RX(c,1) == sum(Para.Interval_len(tempno,:) .* ...
%                     (Var.scenario(c).QXmean(tempno,:) - Para.sign(c) * Var.RelaxIC)) - Var.RelaxIR];
%                 Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Var.RX(c,1)];
%             end
%         else
%             Cons = [Cons, Var.MX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Para.qmax(tempno) * Var.dual_mumax(tempno) ...
%                                 + Para.qmin(tempno) *  Var.dual_mumin(tempno)];
%             Cons = [Cons, Var.RX(c,1) == Para.xpara(tempno) * Var.QX(tempno) - Var.MX(c,1)];      
%         end 
%         
                

        if bidPara.Cstrategic(c) == 1
            if length(tempno) > 1 % a company with multiple units
                temp_len = sum(sum(bidPara.Interval_len(tempno,:) .^2,1).^0.5); %这个是用来计算RelaxIC钱的项目
                if isfield(Setting,'new_version') && Setting.new_version == 1
                    % 如果采用new_version1, 则此时要让Qmax和RelaxIC/RelaxIR挂钩
                    Result.bidRX(c,1) = sum(sum(bidPara.Interval_len(tempno,:) .* ...
                        Result.scenario(c).QXmean(tempno,:))) - Result.kIC * ...
                        temp_len * sum(bidPara.qabsmax(tempno)) - Result.kIR * sum(bidPara.qabsmax(tempno));
                elseif isfield(Setting,'new_version') && Setting.new_version == 2
                    % new_version2这里的temp_len是个向量，针对每个p去计算的
                    % 这里的Para.Interval_len本来就是有正有负，对于买方主体而言是正，对于卖方主体而言是负。    
                    %那么对于卖方主体而言，QXmean与ParaInterval相乘本来就是正数
                    %templen一定是正数
                    temp_len = sum(bidPara.Interval_len(tempno,:) .^2,1).^0.5; 
                    Result.bidRX(c,1) = sum(max(0,sum(bidPara.Interval_len(tempno,:) .* ...
                        Result.scenario(c).QXmean(tempno,:),1) - Result.kIC * sum(bidPara.qabsmax(tempno)) * temp_len)) ... 
                         - Result.kIR * sum(bidPara.qabsmax(tempno));
                    
                elseif isfield(Setting,'new_version') && Setting.new_version == 11
                    Result.bidRX(c,1) = sum(max(0,sum(bidPara.Interval_len(tempno,:) .* ...
                        Result.scenario(c).QXmean(tempno,:),1) - Result.kIC * abs(bidPara.qabsmax(tempno)) * abs(bidPara.Interval_len(tempno,:)) ...
                        ))  - Result.kIR * sum(bidPara.qabsmax(tempno));

                else
                    % 如果采用的不是新version，此时的relaxIC, relaxIR就不会和Qmax结合在一起
                    Result.bidRX(c,1) = sum(sum(bidPara.Interval_len(tempno,:) .* ...
                        Result.scenario(c).QXmean(tempno,:))) - Result.RelaxIC * temp_len - Result.RelaxIR;
                end
%                 Result.MX(c,1) = Para.xpara(tempno) * Result.QX(tempno) - Result.bidRX(c,1);
            else % a single unit
                % 为什么此时要加入Para.sign。因为此时的Para.Interval_len是在括号外相乘
                if isfield(Setting,'new_version') && Setting.new_version == 1
                    % 此时引入了和Qmax的相关项
                    Result.bidRX(c,1) = sum(bidPara.Interval_len(tempno,:) .* ...
                        (Result.scenario(c).QXmean(tempno,:) - bidPara.sign(c) * Result.kIC * sum(bidPara.qabsmax(tempno)))) - Result.kIR * sum(Para.qabsmax(tempno));
                elseif isfield(Setting,'new_version') && (Setting.new_version == 2 || Setting.new_version == 11)
                    % 此时考虑了和max项目
                    Result.bidRX(c,1) = sum(max(0,bidPara.Interval_len(tempno,:) .* ...
                        (Result.scenario(c).QXmean(tempno,:) - bidPara.sign(c) * Result.kIC * ...
                        sum(bidPara.qabsmax(tempno))))) - Result.kIR * sum(bidPara.qabsmax(tempno));
                else
                    % 如果是最原始的模式，relaxIC,relaxIR与Qmax无关。此时如果是买方，Para.sign=1、卖方则反之
                    Result.bidRX(c,1) = sum(bidPara.Interval_len(tempno,:) .* ...
                        (Result.scenario(c).QXmean(tempno,:) - bidPara.sign(c) * Result.RelaxIC)) - Result.RelaxIR;
                end
%                 Result.MX(c,1) = Para.xpara(tempno) * Result.QX(tempno) - Result.bidRX(c,1);
            end
 
            Result.MX(c,1) = bidPara.xpara(tempno) * Result.QX(tempno) - Result.bidRX(c,1);
            Result.RX(c,1) = Para.xpara(tempno) * Result.QX(tempno) - Result.MX(c,1);
            
%             if isfield(Setting,'falsebid_volume') %增加falsebid
%                 % QX表示得到物品，MX表示付出钱，RX表示效用
%                 Result.MX(c,1) = (Para.xpara(tempno) + Setting.falsebid_volume(tempno)) * Result.QX(tempno) - Result.bidRX(c,1);
%                 Result.bidRX(c,1) = Result.bidRX(c,1) - Setting.falsebid_volume(tempno) * Result.QX(tempno); 
%             else
%                 Result.MX(c,1) = bidPara.xpara(tempno) * Result.QX(tempno) - Result.bidRX(c,1);
%                 Result.bidRX(c,1) = Result.bidRX(c,1);
%             end
        else
            if isfield(Setting,'new_version') && (Setting.new_version ~= 0)
                % 如果是new_version，还需要考虑到承担kIR部分的值
                Result.MX(c,1) = bidPara.xpara(tempno) * Result.QX(tempno) - bidPara.qmax(tempno) * Result.dual_mumax(tempno) ...
                                + bidPara.qmin(tempno) *  Result.dual_mumin(tempno) + Result.kIR * sum(bidPara.qabsmax(tempno));
            else
                % 如果不是newversion, 直接按照LMP结算
                Result.MX(c,1) = bidPara.xpara(tempno) * Result.QX(tempno) - bidPara.qmax(tempno) * Result.dual_mumax(tempno) ...
                                + bidPara.qmin(tempno) *  Result.dual_mumin(tempno);
            end
            
            % 如果是买家的话，买家出价 > 边际价格，此时达到qmax, 结算价格要略减，没问题
            Result.bidRX(c,1) = bidPara.xpara(tempno) * Result.QX(tempno) - Result.MX(c,1);
            Result.RX(c,1) = Para.xpara(tempno) * Result.QX(tempno) - Result.MX(c,1);
            
            
%             if isfield(Setting,'falsebid_volume') %增加falsebid
%                 % QX表示得到物品，MX表示付出钱，RX表示效用
%                 Result.MX(c,1) = Setting.falsebid_volume(tempno) * Result.QX(tempno) + Result.MX(c,1);
%                 Result.bidRX(c,1) = Result.bidRX(c,1) - Setting.falsebid_volume(tempno) * Result.QX(tempno); 
%             else
%                 Result.bidRX(c,1) = Result.bidRX(c,1);
%             end
        end 
    end
    
    % Calculate BB
    if Setting.topology == 0
        Result.BB =  sum(Result.MX); %收到的所有的钱的集合
    else
        Result.congestion_revenue = sum(Result.dual_rhomax .* Para.BranchMax) + ...
            sum(Result.dual_rhomin .* Para.BranchMax);
        Result.BB = sum(Result.MX) - Result.congestion_revenue; %congestion_revenue是要付出的额外的钱
    end
    
    Result.RelaxBB = -Result.BB;
    Result.surplus = Result.BB;
    
    % 计算消费者福利和生产者福利
    Result.bid_producer_welfare = sum(Result.bidRX(Para.Gset_company));
    Result.bid_consumer_welfare = sum(Result.bidRX(Para.Dset_company));
    Result.producer_welfare = sum(Result.RX(Para.Gset_company));
    Result.consumer_welfare = sum(Result.RX(Para.Dset_company));
    
    % 计算LMP下的各路福利
    if Setting.direct_solve == 0
         Result.LMPI = (Result.LMP + Para.signI' .* Result.spread(:));
%         if isfield(Setting,'settle_spread') && Setting.settle_spread == 0
%             % spread里分为两部分，一部分是规则里的spread，一部分是谎报的量
%             if isfield(Setting,'falsebid_volume') 
%                 Result.LMPI = Result.LMP + Para.signI(:) .* Result.spread(:) + Setting.falsebid_volume(:);
%             else
%                 % 默认全都是谎报的量
%                 Result.LMPI = Result.LMP; %否则按照spread结算
%             end
%         else
%             % sign的话卖家是-1
%             Result.LMPI = (Result.LMP + Para.signI' .* Result.spread(:));  % 如果没有settle spread，其实默认按照LMP减掉申报的spread
%         end
        Result.LMPsettle = Result.LMPI .* Result.QX;
        Result.LMPsettleC = zeros(Num.C,1);
        for c = 1:Num.C
            Result.LMPsettleC(c) = sum(Result.LMPsettle(Para.company_set(c).no));
        end
        Result.Settle_diff = Result.MX - Result.LMPsettleC;
    end
    if Setting.topology == 1 && Setting.direct_solve == 0
        Result.LMPbalance = sum(Result.LMPsettleC) - Result.congestion_revenue;
    end

    
    
    Result.QXPara = zeros(Num.C,1);
    Result.bidQXPara = zeros(Num.C,1);
    for c = 1:Num.C %这个可以补充进来, 但这里指的是基础状态下的生产成本
        tempno = Para.company_set(c).no;
        Result.QXPara(c) = Para.xpara(tempno) * Result.QX(tempno);
        Result.bidQXPara(c) = bidPara.xpara(tempno) * Result.QX(tempno);
    end
    if Setting.direct_solve == 0
        Result.LMPsettleC_RX = Result.QXPara - Result.LMPsettleC;
        Result.LMPsettleC_bidRX = Result.bidQXPara - Result.LMPsettleC;
    end
    

end