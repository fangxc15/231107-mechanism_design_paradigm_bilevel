% clear 
% script_type = 'xls';
% xlsname = 'Clearing_V1';
% sheetname = '8asymmetry';
% Setting_suffix = '4P';
% 
% Setting.NumP = 4;
% Setting.sample_method = 1; %0表示采用等距离抽样,1表示采用混合抽样,2表示采用等概率增量抽样
% [Para,Num] = read_info([xlsname,'.xlsx'],sheetname, Setting);
% 
% 
% Para.SW_max = 0;%这个离散化了以后会有很多问题
% Setting.Qipp = 1;
% Setting.SWtolerance = 1e-5;
% Setting.bigM = 10000;

%
%%
max_spread = 0.5;
choose_spread = 0.25;


tolerance = 3; %表示几位小数
resolution = 1/10^tolerance;
start_point = round(min(Para.xmin - max_spread),2);
end_point = round(max(Para.xmax + max_spread),2);
temp_point = start_point:resolution:end_point;

% p_cumu代表累计分布曲线
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

%%
Gnum = length(Para.Gset);
Dnum = length(Para.Dset);
totalnum = Gnum+Dnum;
% clear coefs
% % coefs的行是个数，列是价格。Value代表低于该成本价的主体数量的概率
% for p = 1:size(p_cumu,2)
%     coefs(:,p) = cal_coefs(p_cumu(:,p));
% end
%%

%对于需求侧来讲出清的概率
for i = 1:Num.I
    
    temp_Gset = setdiff(Para.Gset,i);
    temp_Dset = setdiff(Para.Dset,i);
    
    % 选出需求侧的系数和累计概率
    if Para.type(i) == 1 %如果是买方
%         temp_coefs = coefs(:, (Para.xmin(i) - start_point - max_spread)/resolution + 1:(Para.xmax(i) - start_point)/resolution + 1);
        temp_cumu = p_cumu(:, (Para.xmin(i) - start_point - max_spread)/resolution + 1:(Para.xmax(i) - start_point)/resolution + 1);
        G_setoff = -choose_spread/resolution;
        D_setoff = 0;
    else %如果是卖方
%         temp_coefs = coefs(:, (Para.xmin(i) - start_point)/resolution + 1:(Para.xmax(i) - start_point + max_spread)/resolution + 1);
        temp_cumu = p_cumu(:, (Para.xmin(i) - start_point)/resolution + 1:(Para.xmax(i) - start_point + max_spread)/resolution + 1);
        G_setoff = 0;
        D_setoff = choose_spread/resolution;
    end
    
    % 从price_spread到最终的累计概率曲线，好像不是一个线性关系，没办法很好的建模.
    % 如果我要是直接拿一个固定的去分析，那起不到优化的效果。是否可以去逼近？
    
    
    
    for p = 1:size(temp_cumu,2)
        temp_prob = temp_cumu(i,p);
        input_vector = temp_cumu([1:i-1 i+1:Num.I],p);
        
        
        collect(i).coefs(:,p) = [cal_coefs(input_vector);0];
%         % 计算每个主体的剩余累计概率
%         if temp_prob < 0.02
%             collect(i).coefs(:,p) = temp_coefs(:,p);
%         elseif 1-temp_prob < 0.02
%             collect(i).coefs(:,p) = [temp_coefs(2:Num.I+1,p);0]; 
%         else
%             diag_matrix = diag((1-temp_prob) * ones(1,length(temp_cumu(:,p))+1),0) + diag(temp_prob *  ones(1,length(temp_cumu(:,p))),-1);
%             if abs(det(diag_matrix)) < 1e-3
%                 fprintf('fail %d %d %f %f\n',i,p,det(diag_matrix),temp_prob);
%             else
%                 fprintf('success %d %d %f %f\n',i,p,det(diag_matrix),temp_prob);
%             end 
%             collect(i).coefs(:,p) = inv(diag_matrix) * temp_coefs(:,p);
%         end
    end
%     collect(i).coefs(Num.I+1,:) = zeros(1,size(temp_cumu,2));
end



%%
% 计算期望出清量：对于Demand而言
% 对于D而言,能得到出清的时候是有5个小于它。因为g>D-1-d, g+d>=D; 否则，不能得到出清
for i = 1:length(Para.Dset)
    ii = Para.Dset(i);
    Exp_q_raw(ii,:) = sum(collect(ii).coefs(Dnum+1:totalnum+1,:));
end

% 计算期望出清量：对于Generator而言
% 对于G而言，得不到出清的状态是：因为g>=D-d, g+d>=D;  得到出清的状态是g<D-d, g+d<D
for i = 1:length(Para.Gset)
    ii = Para.Gset(i);
%     Exp_q(ii,:) = -sum(collect(ii).cumu_p(1:Dnum,:));
    Exp_q_raw(ii,:) = -1+sum(collect(ii).coefs(Dnum+1:totalnum+1,:));
end
%%
Gend_point = (max_spread - choose_spread)/resolution + 1;

for i = 1:Num.I
    if Para.type(i) == 1
        Dstart_point = 1 + max_spread/resolution - choose_spread/resolution;
        Dend_point = size(Exp_q_raw(i,:),2) - choose_spread/resolution;
        Exp_q(i,:) = Exp_q_raw(i,Dstart_point:Dend_point);
    elseif Para.type(i) == 2
        Gstart_point = 1 + choose_spread/resolution;
        Gend_point = size(Exp_q_raw(i,:),2) - max_spread/resolution + choose_spread/resolution;
        Exp_q(i,:) = Exp_q_raw(i,Gstart_point:Gend_point);
    end
end

%%
new_Para = Para;
newNum = Num;
% new_Para.xmax = Para.xmax;
% new_Para.xmin = Para.xmin;
% new_Para.distribution = Para.distribution;
% new_Para.type = Para.type;

% 只有new_Para的point要改
new_Para.Point = [];
for i = 1:newNum.I
    if new_Para.type(i) == 1 %如果是买方
        new_Para.Point(i,:) = new_Para.xmin(i):resolution:new_Para.xmax(i);
    else
        new_Para.Point(i,:) = new_Para.xmin(i):resolution:new_Para.xmax(i);
    end
        
end


newNum.P = (new_Para.xmax(i) - new_Para.xmin(i))/resolution;
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

%%
for i = 1:newNum.I
   
    for p = 1:newNum.P+1
        if Para.type(i) == 1 %如果是买方
            rip(i,p) = sum(sum((Exp_q(i,1:p-1)+Exp_q(i,2:p))/2 .* new_Para.Interval_len(i,1:p-1)));
%              Cons = [Cons, Var.rip(i,p) == sum(sum((Var.qip_mean(i,1:p-1) - Var.RelaxIC) .* Para.Interval_len(i,1:p-1)))];
%              Cons = [Cons, Var.mip(i,p) == Para.Point(i,p) * Var.qip(i,p) + Var.RelaxIR - Var.rip(i,p)];
        end
        if Para.type(i) == 2 %如果是卖方
            rip(i,p) = sum(sum(-(Exp_q(i,p+1:newNum.P+1) + Exp_q(i,p:newNum.P))/2 .* new_Para.Interval_len(i,p:newNum.P)));
%              Cons = [Cons, Var.rip(i,p) == sum(sum((-Var.qip_mean(i,p:Num.P) - Var.RelaxIC) .* Para.Interval_len(i,p:Num.P)))];
%              Cons = [Cons, Var.mip(i,p) == Para.Point(i,p) * Var.qip(i,p) + Var.RelaxIR - Var.rip(i,p)];
        end
    end
end
%% 每个主体的期望出清量
% choose_spread

% 买方和卖方不一样,如果认为是降低0.1？

%%
fprintf('期望交割量')
sum(new_Para.Prob .* Exp_q,2)
fprintf('期望总交割量')
sum(sum(new_Para.Prob .* Exp_q,2))

fprintf('期望社会福利')
sum(sum(new_Para.Point.* new_Para.Prob .* Exp_q))

% 仅在Groves机制类内，看能否收支平衡, 怎么看? 因为出清目标是社会福利最大化, 
% Groves机制类为了能确保激励相容,给他的效用必然是它能创造出的社会福利增量
% 因此我们只要看它的效用即可。同时我们又知道，满足激励相容的效用函数，必然是q(x)的积分
% 因此在社会福利最大化的出清规则下, 为了保证激励相容，q(x)的积分就是效用。同时这个一定是社会福利增量

% 如果本身的机制目标不是社会福利最大化
% 那q(x)的积分仅仅是它的效用而已,并不是社会福利的增量. 而且这个q(x)相当于为了保证激励相容所需要付出的最小的支付量.
 
fprintf('每个主体获得效的效用')
sum(new_Para.Prob .* rip,2)
fprintf('主体获得的效用总额')
sum(sum(new_Para.Prob .* rip,2))

% fprintf('期望社会福利增量')

%%
% for i = 1:10000
%     temp_prob = i/10000;
%     temp_matrix = diag((1-temp_prob) * ones(1,length(p_cumu(:,p))+1),0) + diag(temp_prob *  ones(1,length(p_cumu(:,p))),-1);
%     det_matrix(i) = det(temp_matrix);
% end
%%
% 有没有办法把概率还原回去
% 有4件发生的概率
% 有3件发生的概率

% 排除了某个之后，有4件发生的概率应该是，有5件发生的概率

% d
% y6 = p*x5
% y5 = p*x4 + (1-p)x5
% y4 = p*x3 + (1-p)x4;
% 
% y1 = p*x0 + (1-p)x1;
% y0 =        (1-p)x0;
% 
% x5 = y6/p
% x4 = (y5 - (1-p)x5)/p
% x4 = (y5 - (1-p)y6/p)p
%    = y5/p - (1-p)/p^2*y6
% x3 = (y4 - (1-p)x4)/p
%    = y4/p - y5(1-p)/p^2 - (1-p)^2/p^3 y6
% 
% 
    

% for i1 = 1:Num.I
%     point = linspace(Para.xmin(i1),Para.xmax(i1),Num.P2);
%     for i2 = 1:Num.I
%         Para.
%     end
% end 
%%
% coefs = cal_coefs([0.4;0.7])