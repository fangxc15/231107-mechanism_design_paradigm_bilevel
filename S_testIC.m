%% 这里是测试主体会怎么谎报

Setting.sample_method = 1;
Setting.NumP = 50;
Num.P = 50;
[Para,Num] = process_Para(Para,Num,Setting);

falserate_array = 0:0.05:0.5;
kIC_array = 0:0.02:0.3;

% 这里首先去研究不同的主体、在不同的kIC下，不同策略的收益
for strategic_company = 1:Num.C
%     strategic_company=2;
    input_spread = zeros(1,Num.I);
    Setting.falsebid_volume = zeros(1,Num.I);
    strategic_unit = Para.company_set(strategic_company).no;
    % strategic_volume = 3;
    % target_price = 38.7;
    falserate_no = 0;
    clear Result_collect
    for falserate_no = 1:length(falserate_array)%0:0.05:0.5
        false_rate = falserate_array(falserate_no);
    % false_rate = 0.4;
        target_bidding = max(Para.xpara(strategic_unit)) * (1+false_rate);
        Setting.falsebid_volume(strategic_unit) = target_bidding - Para.xpara(strategic_unit); %Para.xpara(strategic_unit) * false_rate
%         falserate_no = falserate_no + 1;
        kIR = 0;
        kIC_no = 1;
        kIC = kIC_array(kIC_no);
        
    %     Setting.settle_spread = 0;
    %     Setting.falsebid_volume = input_spread;
        Result_falsebid = optimization_primal(Num,Para,Setting, input_spread, kIR, kIC);
        Result_falsebid.kIC = kIC;
        Result_falsebid.false_rate = false_rate;
        Result_falsebid.strategic_company =  strategic_company;
        Result_collect(falserate_no, kIC_no) = Result_falsebid;
    %     Result_collect(falserate_no, kIC_no).kIC = kIC;
    %     Result_collect(falserate_no, kIC_no).false_rate = false_rate;
    %     Result_collect(falserate_no, kIC_no).strategic_company = strategic_company;

        for kIC_no = 2:length(kIC_array)%0.02:0.02:0.3
%             kIC_no = kIC_no + 1;
            kIC = kIC_array(kIC_no);
            Result_falsebid.kIC = kIC;
            Result_falsebid.false_rate = false_rate;
            Result_falsebid.strategic_company =  strategic_company;

            Result_falsebid = calculate_MX(Result_falsebid, Setting, Num, Para);
            Result_falsebid = calculate_Relax(Result_falsebid, Setting, Num, Para);
            Result_collect(falserate_no, kIC_no) = Result_falsebid;
    %         Result_collect(falserate_no, kIC_no).kIC = kIC;
    %         Result_collect(falserate_no, kIC_no).false_rate = false_rate;
    %         Result_collect(falserate_no, kIC_no).strategic_company = strategic_company;
        end
    end
    
    % Result_collect记录了不同的谎报幅度/kIC幅度下市场结果，市场结果汇总形成各个matrix
    
    RX_matrix = zeros(size(Result_collect,1), size(Result_collect,2));
    bidRX_matrix = zeros(size(Result_collect,1), size(Result_collect,2));
    QX_matrix = zeros(size(Result_collect,1), size(Result_collect,2));
    falserate_matrix = zeros(size(Result_collect,1), size(Result_collect,2));
    kIC_matrix = zeros(size(Result_collect,1), size(Result_collect,2));
    LMPsettleC_RX = zeros(size(Result_collect,1), size(Result_collect,2));


    for i = 1:size(Result_collect,1)
        for j = 1:size(Result_collect,2)
            RX_matrix(i,j) = Result_collect(i,j).RX(strategic_company);
            bidRX_matrix(i,j) = Result_collect(i,j).bidRX(strategic_company);
            QX_matrix(i,j) = sum(Result_collect(i,j).QX(strategic_unit));
            falserate_matrix(i,j) = Result_collect(i,j).false_rate;
            kIC_matrix(i,j) = Result_collect(i,j).kIC;
            LMPsettleC_RX(i,j) = Result_collect(i,j).LMPsettleC_RX(strategic_company);
        end
    end

    tempcollect.RX_matrix = RX_matrix;
    tempcollect.LMPsettleC_RX = LMPsettleC_RX;
    tempcollect.QX_matrix = QX_matrix;
    tempcollect.bidRX_matrix = bidRX_matrix;

    Result_parti_collect(strategic_company) = tempcollect;
end

% 各个matrix汇总放在result_parti_collect里，表征各主体在不同的kIC之下会谎报多少
%% 这里是为了看主体的最优策略是什么，放在Plot_index里
index = zeros(Num.C,size(RX_matrix,2));
for c = 1:Num.C  %8/10不是策略性的
    tempno = Para.company_set(c).no;
    Result_parti_collect(c).RX_ratio_matrix = (Result_parti_collect(c).RX_matrix - Result_parti_collect(c).RX_matrix(1,:))/sum(Para.xpara(tempno).* Para.qabsmax(tempno));
    % Result_parti_collect(c).RX_matrix./Result_parti_collect(c).RX_matrix(1,:);
    RX_ratio_matrix = Result_parti_collect(c).RX_ratio_matrix;
    a = size(RX_ratio_matrix,1);
    b = size(RX_ratio_matrix,2);
    RX_ratio_matrix(2:a,:) = RX_ratio_matrix(2:a,:) - 0.05;
    [~, index(c,:)] = max(RX_ratio_matrix);
end

plot_index = falserate_array(index);
%% 这里是为了画出在不同的kIC下主体的最优谎报幅度
set(groot,'defaultfigurePosition',[200 200 340 280]);
%     set(groot,'defaultLegendFontName','Times New Roman');
set(groot,'defaultLegendFontSize',12);
set(groot,'defaultAxesFontSize',11);
%     set(groot,'defaultFontSize',14);

set(groot,'defaultAxesFontWeight','bold');
set(groot,'defaultAxesFontName','Times New Roman');
set(groot,'defaultAxesFontName',['SimSun']);
set(0,'defaultfigurecolor','w'); %设置背景颜色为白色
% xvalues = {0:0.02:0.3};
xvalues = {};
for kIC_no = 1:length(kIC_array)
    xvalues(kIC_no) = {num2str(kIC_array(kIC_no),'%.2f')} ;
end
% xvalues = num2cell(0:0.02:0.3);
yvalues = {'卖家1','卖家2','卖家3','卖家4','买家5','买家6','买家7','买家8','买家9','买家10'};
figure(1)
h = heatmap(xvalues,yvalues,plot_index);
h.FontName = ['SimSun'];
h.HandleVisibility = 'on';
h.CellLabelColor = 'none';
% h.CellLabelFormat = ''
% set(gca,'FontWeight','bold')
xlabel('激励相容松弛度')
% set(gca,'xticklabels',{0:0.02:0.3})
ylabel('主体编号')
title('虚报幅度')
ifsave = 1;
if ifsave
    print('-dpng','-r1000',[Picture_folder,'/','participant_strategy.png']);
    saveas(1,[Picture_folder,'/','participant_strategy.jpg'])
end

%% 可以查看在不同的kIC之下实现的社会福利是多少？
input_spread = zeros(1,Num.I);
kIR = 0;
for kIC_no = 1:length(kIC_array)
    kIC = kIC_array(kIC_no);
    Setting.falsebid_volume = zeros(1,Num.I);
    for strategic_company = 1:Num.C
        strategic_unit = Para.company_set(strategic_company).no;
        if plot_index(strategic_company,kIC_no) > 0 
            target_bidding = max(Para.xpara(strategic_unit)) * (plot_index(strategic_company,kIC_no) + 1);
            Setting.falsebid_volume(strategic_unit) = target_bidding - Para.xpara(strategic_unit); %Para.xpara(strategic_unit) * false_rate
        end
    end
    Result_bid = optimization_primal(Num,Para,Setting, input_spread, 0, kIC);
    Result_kIC_collect(kIC_no) = Result_bid;
end
[Result_kIC_collect.SW] %这个相当于是展现不同的kIC幅度下实现的社会福利差异

mkdir('Result_test_false')
save(['Result_test_false/test_data.mat'])
%% 所以这个到底要实现一个什么功能呢？目的感觉很不明确。
% 这个感觉没必要写的这么复杂，因为到此为止就写这些感觉足够了。
% 是否要根据它的需求选择一个点？

% 可以评估各个分立主体的上限

%% 
% 一个思路是评估各个主体在不同的kIC下，虚假申报能获取多少超额，如果没有就说明这个kIC可被接受
% 另一个思路是要选取一个给定的合理的值