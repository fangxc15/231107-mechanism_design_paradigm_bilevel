function [Para,Num] = process_Para(Para,Num,Setting)
    
    % 可以提出一个抽样算法, 可能会比较复杂
    % 抽样主要考虑等距离增加即可(在决策spread的时候).如果考虑决策Q(X),也可以考虑等概率密度增加.
        % 而如果等距离, 会有利于计算M(X). 另外, 如果在决策spread时候，没有概率的概念
        % 如果等概率密度增加，或许在Q(X)的时候计算
            %等概率密度增加，在x_better一侧会有一些问题，这个时候必然会选社会福利最大的。对这种情况怎么处理呢？
        % 而距离过大不能过大，这会影响用Q(x_mid)代表Q(x)的准确性,同时影响r(x)[因为要对q(x)做积分]项的准确性]
    % sample_method:0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样

    Para.xmax = Para.xmax(:)';
    Para.xmin = Para.xmin(:)';
    Para.type = Para.type(:)';
    Para.qmax = Para.qmax(:)';
    Para.qmin = Para.qmin(:)';
    Para.company = Para.company(:)';
    Para.distribution = Para.distribution(:)';
    

%     %%
% %     % 这里的Normpoint是交替上升的. 其实新的抽样方式, 应当是根据Para.xpara来确定的.
% %     Para.normPoint = ones(Num.I,1) * [0:Num.I:Num.I*Num.P]/(Num.I*(Num.P+1)-1) + [0:Num.I-1]'/(Num.I*(Num.P+1)-1);
%     Para.Point = (Para.xmax - Para.xmin)' * ones(1,Num.P+1) .* Para.normPoint + Para.xmin';
% 
%     if Para.distribution(i) == 1
%         norm_mid = (Para.xmax + Para.xmin)/2;
%         norm_std = (Para.xmax - Para.xmin)/6;
%         point_temp = repmat(norm_std',1,Num.P+1) .* norminv(Para.normPoint * (normcdf(3)-normcdf(-3)) + normcdf(-3)) + norm_mid';
%         if Setting.sample_method == 1
%             Para.Point = (Para.Point + point_temp)/2;
%         elseif Setting.sample_method == 2
%             Para.Point = point_temp;
%         end
%     end
%     
%     for i = 1:Num.I
%         % 如果是均匀分布, 等距离/等概率密度增加是一致的.
% %         Para.Point(i,:) = linspace(Para.xmin(i),Para.xmax(i),Num.P+1); %:Para.xmax(i)
%         if Para.distribution(i) == 1
%             norm_mid = (Para.xmax(i) + Para.xmin(i))/2;
%             norm_std = (Para.xmax(i) - Para.xmin(i))/6;
%             point_temp = norm_std * norminv(Para.normPoint(i,:) * (normcdf(3)-normcdf(-3)) + normcdf(-3)) + norm_mid;
% %             point_temp = (norminv(linspace(normcdf(-3),normcdf(3),Num.P+1))) * norm_std + norm_mid;
%             if Setting.sample_method == 1     %两种抽样方式混合
%                 Para.Point(i,:) = (Para.Point(i,:) + point_temp)/2;
%             elseif Setting.sample_method == 2 %纯按照概率抽样
%                 Para.Point(i,:) = point_temp;
%             end
%         end 
%     end 
% 
%     Para.Interval_len = [Para.Point(:,1)-Para.xmin' Para.Point(:,2:Num.P+1)-Para.Point(:,1:Num.P) Para.xmax'-Para.Point(:,Num.P+1)];
%     这个Probability假设的是均匀分布. 如果从r(x)的积分角度来看, 应该采用等距离积分. 到底采用啥积分，不知道.
%     for i = 1:Num.I
%         if Para.distribution(i) == 0 
%             Para.p_cumu(i,:) = (Para.Point(i,:) - Para.xmin(i))/(Para.xmax(i) - Para.xmin(i));
%         elseif Para.distribution(i) == 1
%             norm_mid = (Para.xmax(i) + Para.xmin(i))/2;
%             norm_std = (Para.xmax(i) - Para.xmin(i))/6;
%             Para.p_cumu(i,:) = normcdf((Para.Point(i,:) - norm_mid)/norm_std);
%         end
%     end
% 
%     Para.Prob = 0.5 * ([Para.p_cumu(:,2:Num.P+1) 2-Para.p_cumu(:,Num.P+1)] - [-Para.p_cumu(:,1) Para.p_cumu(:,1:Num.P)]); 
    
    Para.Dset = find(Para.type == 1); %1是买方
    Para.Gset = find(Para.type == 2); %2是卖方
    
    Para.signI = zeros(1,Num.I);
    Para.signI(Para.Dset) = 1; %买方1
    Para.signI(Para.Gset) = -1; %卖方-1

    
    Para.totalD = sum(Para.qmax(find(Para.type == 1)));
    Para.totalG = -sum(Para.qmin(find(Para.type == 2)));
    for c = 1:Num.C
        tempno = find(Para.company == c);
        Para.company_set(c).no = tempno;
        if Para.type(tempno(1)) == 1 
           Para.sign(c) = 1; %买方1
           if isfield(Setting,'percent_strategic')
               if sum(Para.qmax(tempno)) > Para.totalD * Setting.percent_strategic
                   Para.Cstrategic(c) = 1;
               else
                   Para.Cstrategic(c) = 0;
               end
           else
               Para.Cstrategic(c) = 1;
           end
        else 
           Para.sign(c) = -1;%卖方-1
           if isfield(Setting,'percent_strategic')
               if -sum(Para.qmin(tempno)) > Para.totalG * 0.1
                   Para.Cstrategic(c) = 1;
               else
                   Para.Cstrategic(c) = 0;
               end
           else
               Para.Cstrategic(c) = 1;
           end
        end   
        %Para.SIC_sign = Para.type * 2 - 3;%假设是买方则为负号，假设为卖方则为正号
    end
    
    Para.Dset_company = find(Para.sign == 1);
    Para.Gset_company = find(Para.sign == -1);
%     if isfield(Setting,'percent_strategic')
%         strategic_Cno = find(Para.Cstrategic == 1);
%     else
%         strategic_Cno = 1:Num.C;
%     end
    
    % 这个我希望是只有strategic的部分才有probability
%     Para.prob = tempprob_c/(sum(sum(tempprob_c)) - sum(tempprob_c(2:Num.C,1)));% 这个是计算SW时候的概率信息, c行p列
    if isfield(Setting,'sample_method') && Setting.sample_method == 2
       Para.normPoint = zeros(Num.I,Num.P+1);
        for i = 1:Num.I
            if Para.type(i) == 1
                Para.normPoint(i,:) = [0.5,0.5-0.5/(Num.P+1)*(1.5:Num.P+0.5)]; 
            elseif Para.type(i) == 2
                Para.normPoint(i,:) = [0.5,0.5+0.5/(Num.P+1)*(1.5:Num.P+0.5)]; 
            end
        end
        for i = 1:Num.I
            Para.Point(i,:) = norminv(Para.normPoint(i,:)) * Para.xstd(i) + Para.xpara(i);
        end
        tempprob_i = 1/(Num.P+1) * ones(Num.I,Num.P+1);
    else
        % 这个地方为什么可以这么做？从正常逐渐逼近worse
        Para.Point = -(Para.xpara - Para.xworse)' * [0:Num.P]/Num.P + Para.xpara'; %这个是仅限spread决策法的时候;
        tempprob_i = normpdf((Para.Point - Para.Point(:,1)) ./repmat(Para.xstd',1,Num.P+1));
    end
    Para.Interval_len = Para.Point(:,1:Num.P) - Para.Point(:,2:Num.P+1);

        
%     Para.mid_Point = (Para.Point(:,1:Num.P-1 ) + Para.Point(:,2:Num.P))/2;
%     for i = 1:Num.I
%         Para.mid_norm(i,:) = (Para.mid_Point(i,:) - Para.Point(i,1))/Para.xstd(i);
%         tempprob_i(i,1) = abs(normcdf(Para.mid_Point(i,1)) - 0.5) * 2;
%         for p = 2:Num.P
%             tempprob_i(i,p) = abs(normcdf(Para.mid_Point(i,p)) - normcdf(Para.mid_Point(i,p-1)));
%         end
%         tempprob_i
%     end
%     Para.prob
%     


    
    
    
    for c = 1:Num.C
        
        if Para.Cstrategic(c) == 1
            tempprob_c(c,:) = mean(tempprob_i(Para.company_set(c).no, :),1);
        else
            tempprob_c(c,:) = zeros(1,Num.P+1);
        end
    end
    Para.prob = tempprob_c / (sum(sum(tempprob_c(:,2:Num.P+1))) + max(tempprob_c(:,1)));
    
    Para.qabsmax = zeros(1,Num.I);
    Para.qabsmax(Para.Dset) = Para.qmax(Para.Dset);
    Para.qabsmax(Para.Gset) = abs(Para.qmin(Para.Gset));
    for c = 1:Num.C
        Para.qabsmax_C(c) = sum(Para.qabsmax(Para.company_set(c).no));
    end
    
end