function Result = calculate_Relax(Result, Setting, Num, Para)

    % 这里应该也分为真实的和虚假的
    
    if Result.unified_spread == 0
        Para.SW_max = Result.SW;
        Result.RelaxSW = 0;
    else
        Result.RelaxSW = Para.SW_max - Result.SW;
    end

    Result.RelaxIC = Result.kIC/min(Para.xpara);
    Result.RelaxIC = max(0,Result.RelaxIC);

    for c = 1:Num.C
        tempno = Para.company_set(c).no;
        Result.RelaxIR_C(c) = max(0,-Result.RX(c)/ sum(Para.qabsmax(tempno).* Para.xpara(tempno)));
    end
    Result.RelaxIR = sum(Result.RelaxIR_C);
end