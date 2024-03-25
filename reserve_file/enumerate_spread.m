function [result_collect_spread,new_Para] = enumerate_spread(Para,Num,choose_spread_vector,tolerance)
    % 


    %choose_spread_vector: 穷举了vector是多少
    %tolerance: 精度的忍耐度
    max_spread = max(abs(choose_spread_vector));
    resolution = 1/10^tolerance;
    
    %% 构建new_Para, new_Para里的Point是每个主体私有的Point
    new_Para = Para;
    newNum = Num;

    new_Para.Point = [];
    for i = 1:newNum.I
            new_Para.Point(i,:) = new_Para.xmin(i):resolution:new_Para.xmax(i);
    end

    % 这里的new Para
    newNum.P = (new_Para.xmax(i) - new_Para.xmin(i))/resolution;
    %这里的newNum本来就挺有问题的, 每个主体可能算出来不太一样
    
    % 这个地方的new_Para的Interval_len没必要改,因为这里的头尾刚好落在边界上了已经.
    new_Para.Interval_len = new_Para.Point(:,2:newNum.P+1) - new_Para.Point(:,1:newNum.P);
    % 这个Probability假设的是均匀分布. 如果从r(x)的积分角度来看, 应该采用等距离积分. 到底采用啥积分，不知道.
    new_Para.p_cumu = [];
    for i = 1:newNum.I
        if new_Para.distribution(i) == 0 
            new_Para.p_cumu(i,:) = (new_Para.Point(i,:) - new_Para.xmin(i))/(new_Para.xmax(i) - new_Para.xmin(i));
        elseif new_Para.distribution(i) == 1
            norm_mid = (new_Para.xmax(i) + new_Para.xmin(i))/2;
            norm_std = (new_Para.xmax(i) - new_Para.xmin(i))/6;
            new_Para.p_cumu(i,:) = normcdf((new_Para.Point(i,:) - norm_mid)/norm_std);
        end
    end
    new_Para.p_cumu = min(max(0,new_Para.p_cumu),1);
    new_Para.Prob = 0.5 * ([new_Para.p_cumu(:,2:newNum.P+1) 2-new_Para.p_cumu(:,newNum.P+1)] - [-new_Para.p_cumu(:,1) new_Para.p_cumu(:,1:newNum.P)]); 

        
    %% 做spread_num的循环
    for spread_num = 1:length(choose_spread_vector)
        % temp_point是公用的，并且通过一个近似的估计和主体的Para建立联系
        spread_num
        choose_spread = choose_spread_vector(spread_num);
        start_point = round(min(Para.xmin - abs(choose_spread)),2);
        end_point = round(max(Para.xmax + abs(choose_spread)),2);
        temp_point = start_point:resolution:end_point;

        % p_cumu代表累计分布曲线
        clear p_cumu
        clear collect
        clear Exp_q
        
        % 因为计算自己的剩余需求曲线的时候也要用到别人的prob_distribution
        for i = 1:Num.I
            if Para.distribution(i) == 0 
                p_cumu(i,:) = min(max(0,(temp_point - Para.xmin(i))/(Para.xmax(i) - Para.xmin(i))),1);
            elseif Para.distribution(i) == 1
                norm_mid = (Para.xmax(i) + Para.xmin(i))/2;
                norm_std = (Para.xmax(i) - Para.xmin(i))/6;
                p_cumu(i,:) = normcdf((temp_point - norm_mid)/norm_std);
            end
        end

        %p_cumu实际上涵盖了两头延长的部分.
        Gnum = length(Para.Gset);
        Dnum = length(Para.Dset);
        totalnum = Gnum+Dnum;
        % clear coefs
        % % coefs的行是个数，列是价格。Value代表低于该成本价的主体数量的概率
        % for p = 1:size(p_cumu,2)
        %     coefs(:,p) = cal_coefs(p_cumu(:,p));
        % end

        %对于需求侧来讲出清的概率
        %对每个主体依次做一次循环
        for i = 1:Num.I
            
            temp_Gset = setdiff(Para.Gset,i);
            temp_Dset = setdiff(Para.Dset,i);
            % 累计曲线上是通过setoff的设置来考虑price spread
            if Para.type(i) == 1 %如果是买方
                G_setoff = round(-choose_spread/resolution);
                D_setoff = 0;
                self_setoff = round(-choose_spread/resolution);
            else %如果是卖方
                G_setoff = 0;
                D_setoff = round(choose_spread/resolution);
                self_setoff = round(choose_spread/resolution);
            end

            % 从price_spread到最终的累计概率曲线，好像不是一个线性关系，没办法很好的建模成优化.如果我要引入price
            % spread作为变量的话，会产生非线性项
            % 如果我要是直接拿一个固定的去分析，那起不到优化的效果。是否可以去逼近？这里就是穷举逼近的的思路。
            for point = Para.xmin(i):resolution:Para.xmax(i)
                absolute_p = round((point - start_point)/resolution + 1,0);
                relative_p = round((point - Para.xmin(i))/resolution + 1,0);
        %         input_vector = p_cumu([1:i-1 i+1:Num.I],absolute_p + self_setoff);
                %根据这个来求等效(考虑setoff后)的剩余需求曲线(去计算剩余的需求的概率性曲线是多少)
                input_vector = [p_cumu(temp_Gset,absolute_p+G_setoff);p_cumu(temp_Dset,absolute_p+D_setoff)];
                %最后得到的这个是一个概率性曲线
                collect(i).coefs(:,relative_p) = [cal_coefs(input_vector);0];
            end
        end

        % 计算期望出清量：对于Demand而言
        % 对于D而言,能得到出清的时候是有5个小于它。因为g>D-1-d, g+d>=D; 否则，不能得到出清
        for i = 1:length(Para.Dset)
            ii = Para.Dset(i);
            Exp_q(ii,:) = sum(collect(ii).coefs(Dnum+1:totalnum+1,:));
        end

        % 计算期望出清量：对于Generator而言
        % 对于G而言，得不到出清的状态是：因为g>=D-d, g+d>=D;  得到出清的状态是g<D-d, g+d<D
        for i = 1:length(Para.Gset)
            ii = Para.Gset(i);
            Exp_q(ii,:) = -1+sum(collect(ii).coefs(Dnum+1:totalnum+1,:));
        end

        for i = 1:newNum.I
            for p = 1:newNum.P+1
                if new_Para.type(i) == 1 %如果是买方
                    rip(i,p) = sum(sum((Exp_q(i,1:p-1)+Exp_q(i,2:p))/2 .* new_Para.Interval_len(i,1:p-1)));
                end
                if new_Para.type(i) == 2 %如果是卖方
                    rip(i,p) = sum(sum(-(Exp_q(i,p+1:newNum.P+1) + Exp_q(i,p:newNum.P))/2 .* new_Para.Interval_len(i,p:newNum.P)));
                end
            end
        end
        %% 每个主体的期望出清量
        result_collect_spread(spread_num).Exp_q = Exp_q;
        result_collect_spread(spread_num).rip = rip;
        result_collect_spread(spread_num).Exp_q_i = sum(new_Para.Prob .* Exp_q,2);
        result_collect_spread(spread_num).Exp_q_balance = sum(sum(new_Para.Prob .* Exp_q,2));
        result_collect_spread(spread_num).Exp_q_all = sum(sum(max(0,new_Para.Prob .* Exp_q),2));

        result_collect_spread(spread_num).welfare = sum(sum(new_Para.Point.* new_Para.Prob .* Exp_q));
        result_collect_spread(spread_num).welfare_i = sum(new_Para.Prob .* rip,2);
        
        % welfare和分配给它的welfare
        result_collect_spread(spread_num).balance = result_collect_spread(spread_num).welfare - sum(result_collect_spread(spread_num).welfare_i);

        result_collect_spread(spread_num).choose_spread = choose_spread;
        %%
      % 满足激励相容的必然是q(x)的积分. 在Groves机制类内,q(x)的积分就是社会福利的增益.
      % 所以如果是按照社会福利最大化出清的时候, 就看他创造的增量社会福利之和和社会总福利那个高。
      
      % 否则的话，就看分配给他的r(x)之和 与 社会总福利谁高.
      
    %     % 仅在Groves机制类内，看能否收支平衡, 怎么看? 因为出清目标是社会福利最大化, 
    %     % Groves机制类为了能确保激励相容,给他的效用必然是它能创造出的社会福利增量
    %     % 因此我们只要看它的效用即可。同时我们又知道，满足激励相容的效用函数，必然是q(x)的积分
    %     % 因此在社会福利最大化的出清规则下, 为了保证激励相容，q(x)的积分就是效用。同时这个一定是社会福利增量
    % 
    %     % 如果本身的机制目标不是社会福利最大化
    %     % 那q(x)的积分仅仅是它的效用而已,并不是社会福利的增量.
    %     而且这个q(x)相当于为了保证激励相容所需要付出的最小的支付量.[因为在x_worse的时候，会让utilty = 0]

    end
end


