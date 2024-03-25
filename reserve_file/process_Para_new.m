function [Para,Num] = process_Para_new(Para,Num,Setting)
    % 
    % 可以提出一个抽样算法, 可能会比较复杂
    % 抽样既需要考虑等距离, 也需要考虑等概率密度增加, 如何把这两个融在一起?
        % 感觉还是要以等概率密度增加为主，因为会同时影响总的社会福利项、收支平衡项. 给人的感觉还是前者更加重要.
        % 而距离过大，也会影响用x_mid代表x的准确性,同时影响r(x)[因为要对q(x)做积分]项的准确性]
    % sample_method:0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
    %sample_method:0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
    Para.xmax = Para.xmax(:)';
    Para.xmin = Para.xmin(:)';
    Para.type = Para.type(:)';
    Para.qmax = Para.qmax(:)';
    Para.qmin = Para.qmin(:)';
    Para.distribution = Para.distribution(:)';

    Para.normPoint = ones(Num.I,1) * [0:Num.I:Num.I*Num.P]/(Num.I*(Num.P+1)-1) + [0:Num.I-1]'/(Num.I*(Num.P+1)-1);
    Para.Point = (Para.xmax - Para.xmin)' * ones(1,Num.P+1) .* Para.normPoint + Para.xmin';

    for i = 1:Num.I
        % 如果是均匀分布, 等距离/等概率密度增加是一致的.
%         Para.Point(i,:) = linspace(Para.xmin(i),Para.xmax(i),Num.P+1); %:Para.xmax(i)
        if Para.distribution(i) == 1
            norm_mid = (Para.xmax(i) + Para.xmin(i))/2;
            norm_std = (Para.xmax(i) - Para.xmin(i))/6;
            point_temp = norm_std * norminv(Para.normPoint(i,:) * (normcdf(3)-normcdf(-3)) + normcdf(-3)) + norm_mid;
%             point_temp = (norminv(linspace(normcdf(-3),normcdf(3),Num.P+1))) * norm_std + norm_mid;
            if Setting.sample_method == 1     %两种抽样方式混合
                Para.Point(i,:) = (Para.Point(i,:) + point_temp)/2;
            elseif Setting.sample_method == 2 %纯按照概率抽样
                Para.Point(i,:) = point_temp;
            end
        end 
    end 


    Para.Interval_len = Para.Point(:,2:Num.P+1) - Para.Point(:,1:Num.P);
    % 这个Probability假设的是均匀分布. 如果从r(x)的积分角度来看, 应该采用等距离积分. 到底采用啥积分，不知道.
    for i = 1:Num.I
        if Para.distribution(i) == 0 
            Para.p_cumu(i,:) = (Para.Point(i,:) - Para.xmin(i))/(Para.xmax(i) - Para.xmin(i));
        elseif Para.distribution(i) == 1
            norm_mid = (Para.xmax(i) + Para.xmin(i))/2;
            norm_std = (Para.xmax(i) - Para.xmin(i))/6;
            Para.p_cumu(i,:) = normcdf((Para.Point(i,:) - norm_mid)/norm_std);
        end
    end

    Para.Prob = 0.5 * ([Para.p_cumu(:,2:Num.P+1) 2-Para.p_cumu(:,Num.P+1)] - [-Para.p_cumu(:,1) Para.p_cumu(:,1:Num.P)]); 

    Para.Dset = find(Para.type == 1); %1是买方
    Para.Gset = find(Para.type == 2); %2是卖方
    
    if Setting.cumulative_curve == 1
%         resolution = 1/10^(Setting.tolerance);
        for j = 1:Num.I
            for i = 1:Num.I
                if Para.distribution(i) == 0 
                    Para.collect_cumu(j).p_cumu(i,:) = min(max(0,(Para.Point(j,:) - Para.xmin(i))/(Para.xmax(i) - Para.xmin(i))),1);
                elseif Para.distribution(i) == 1
                    norm_mid = (Para.xmax(i) + Para.xmin(i))/2;
                    norm_std = (Para.xmax(i) - Para.xmin(i))/6;
                    Para.collect_cumu(j).p_cumu(i,:) = normcdf((Para.Point(j,:) - norm_mid)/norm_std);
                end
            end
        end
    end
    
end