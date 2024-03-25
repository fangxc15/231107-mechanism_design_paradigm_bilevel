function Result_primal_BBmax = heuristic_find_BBmax(start_spread,end_spread,spread_tolerance,Num,Para,Setting)

    mid_spread = (start_spread + end_spread)/2;
    Result_primal_BBmax(1) = optimization_primal(Num,Para,Setting, mid_spread * ones(1,Num.I));
    i = 1;
    while end_spread - mid_spread > spread_tolerance
        left_mid = (start_spread + mid_spread)/2;
        right_mid = (end_spread + mid_spread)/2;
        i = i + 1;
        Result_primal_BBmax(i) = optimization_primal(Num,Para,Setting, left_mid * ones(1,Num.I));
        i = i + 1;
        Result_primal_BBmax(i) = optimization_primal(Num,Para,Setting, right_mid * ones(1,Num.I));
    %     if record_BB == Result_primal_BBmax(i).BB && (record_BB == Result_primal_BBmax(i-1).BB)
        if Result_primal_BBmax(i).BB > Result_primal_BBmax(i-1).BB %right比较大
            start_spread = mid_spread;
            mid_spread = right_mid;
    %         record_BB = Result_primal_BBmax(i).BB;
        else
            end_spread = mid_spread;
            mid_spread = left_mid;
    %         record_BB = Result_primal_BBmax(i).BB;
        end
    end
end
