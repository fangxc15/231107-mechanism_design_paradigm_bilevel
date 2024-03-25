function Result_primal = heuristic_find_BBbalance(start_spread,end_spread,spread_tolerance,Num,Para,Setting)

    % start_spread = 0;
    % end_spread = 10;
    % spread_tolerance = 0.01;
    clear Result_primal
    Result_primal(1) = optimization_primal(Num,Para,Setting, start_spread * zeros(1,Num.I));
    Result_primal(2) = optimization_primal(Num,Para,Setting, end_spread * ones(1,Num.I));
    start_no = 1;
    end_no = 2;
    mid_spread = (start_spread + end_spread)/2;

    i = 2;
    while end_spread - mid_spread > spread_tolerance
        % 二分搜索的方式,想要搜索BB刚好平衡且SW最大的点。
        % 如果BB>0,说明spread可以往小了搜索一点
        % 如果BB<0,说明spread可以往大了搜索一点? (但这个不一定对, 有可能越往大BB越低)
        % 所以需要设计一个启发式的寻找初值的方法
        input_spread = mid_spread * ones(1,Num.I);
    %     input_spread =  -Para.sign * mid_spread;
        i = i+1;
        Result_primal(i) = optimization_primal(Num,Para,Setting, input_spread);
        if Result_primal(i).BB > 0 % 如果BB>0,说明spread可以往小了搜索一点
            end_spread = mid_spread;
            end_no = i;
            mid_spread = (start_spread + end_spread)/2;
        else % 如果BB<0,说明spread可以往大了搜索一点? (但这个不一定对, 有可能越往大BB越低)
            if Result_primal(start_no).BB > Result_primal(i).BB
                end_spread = mid_spread;
                end_no = i;
                mid_spread = (start_spread + end_spread)/2;

            else %if Result_primal(end_no).BB > Result_primal(i).BB
                start_spread = mid_spread;
                start_no = i;
                mid_spread = (start_spread + end_spread)/2;
            end
        end
    end 
end
%%

