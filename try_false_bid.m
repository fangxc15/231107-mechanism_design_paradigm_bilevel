% 这里一般是直接load某个结果
% 例如load('Result_newversion2//xls_Clearing_IEEE30_Version2_5P_topology1.mat)
%在某个特殊场景下，考虑将某个主体谎报

% 这里spread是机制参数
% Setting.falsebid.volume是它真正谎报的值。

% 这里只是为了考虑让谁谎报
sortrows([(1:Num.I)' Para.signI' Para.xpara' Para.qabsmax' Para.company' Result_maxSW.LMPI Result_maxSW.QX Para.node'],2)
Result_maxSW.QX ./ Para.qabsmax'
Result_falsebid.LMPI
Result_falsebid.QX ./ Para.qabsmax'
%% 为什么不设计成9号主体(买方)具备策略性？因为它谎报造成的社会福利损失不够大
% input_spread = zeros(1,Num.I);
% Setting.falsebid_volume = zeros(1,Num.I);
% strategic_company = 9;
% target_bidding = 35;
% Setting.falsebid_volume([11,12,13]) = [35.8-target_bidding 41.2492-target_bidding 36.9-target_bidding];
%% 2号主体(卖方)距离边际很远，它谎报造成的社会福利损失才是真的很大 0.0978的RelaxIC
% 2号主体的两个机组分别在2号/22号节点，成本28.17/23.49，容量60.97/21.59
% 节点电价分别为35.78/29.35，它全额出清，60.97/21.59
% 可以升高多少而不有损？
% 现在在谎报到38.7后，出清量54.76/7.14，节点电价都是38.7
% %% 设计成2号company

% input_spread = zeros(1,Num.I);
% Setting.falsebid_volume = zeros(1,Num.I);
% strategic_company = 2;
% target_price = 38.7;
% Setting.falsebid_volume([2,3]) = [target_price-28.17 target_price-23.49];

%% 3号主体谎报, 可以算出来0.46/0.91的RelaxIC,边际上的机组本来就是IC非常接近1
input_spread = zeros(1,Num.I);
Setting.falsebid_volume = zeros(1,Num.I);
strategic_company = 3;
target_price = 35;
Setting.falsebid_volume([4]) = [target_price-34.74];



%%
% Setting.settle_spread = 0;
% Setting.falsebid_volume = input_spread;
Result_falsebid = optimization_primal(Num,Para,Setting, input_spread);
% Setting.settle_spread = 1;

%% LMP下真实申报/谎报/VCG机制下的利润
% 对比真实bid和false_bid的情况 %685.31/590.42
Result_falsebid.LMPsettleC_RX(strategic_company) % 这里是按照边际价格结算
Result_maxSW.LMPsettleC_RX(strategic_company) %这里是按照边际价格结算

Result_falsebid.LMPsettleC_bidRX(strategic_company) % 这里是按照边际价格结算
% Result_maxSW.LMPsettleC_bidRX(strategic_company) %这里是按照边际价格结算
Result_maxSW.LMPsettleC_RX(strategic_company) %这里是按照边际价格结算

% VCG机制能获得的利润
Result_vertex(3).RX %和VCG机制下的情况好像不是很一样
Result_vertex(3).RX(strategic_company) %VCG机制下能获取的最大收益是894


%% 社会福利的变化
Result_falsebid.SW  %1391.4, 这个SW好像是视在的SW
Result_maxSW.SW     %1666.1, 社会福利出现大幅下降

%% 社会福利分配（消费者福利？生产者福利的变化）
sum(Result_falsebid.LMPsettleC_RX(Para.Gset_company))  %1134 
sum(Result_maxSW.LMPsettleC_RX(Para.Gset_company))     %851.0

sum(Result_falsebid.LMPsettleC_RX(Para.Dset_company))  %369
sum(Result_maxSW.LMPsettleC_RX(Para.Dset_company))     %613.4

%% baseSW和SW还是有区别的 1697.6/1854.3
Result_falsebid.baseSW
Result_maxSW.baseSW


Result_falsebid.congestion_revenue %193.87
Result_maxSW.congestion_revenue    %389.90

% 851.0+613.4+389.9 = 1854.3
% 1134.0+369+193.9  = 1696.9

% 这里的社会福利计算是对的。因为是通过施加spread来改变出清量。只不过施加spread，此时可能会影响结算价格，但是该影响已经排除。

%% 那LMP机制的RelaxIC到底是多少？
extra_revenue = Result_falsebid.LMPsettleC_RX(strategic_company) - Result_maxSW.LMPsettleC_RX(strategic_company);
Para.qabsmax_C(strategic_company)
norm(Setting.falsebid_volume,2)
extra_revenue/(Para.qabsmax_C(strategic_company)*norm(Setting.falsebid_volume,2)) %0.0621
extra_revenue/(Para.qabsmax * Setting.falsebid_volume')                           %0.0978

Result_vertex(5).RelaxIC %0.1456, 这个是核定的松弛量. 我一定要让我的松弛量比它低就可以。