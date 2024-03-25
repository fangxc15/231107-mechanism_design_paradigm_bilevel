function coefs = cal_coefs(p_cumu)
% 计算剩余供给曲线; 价格越高剩余供给越高. 一次性只能计算一个价格的
% 一开始, 是[1;0;0;0;0;0] 因为一个以上的的概率是0;
% 然后针对每一个新加入的主体，都做一次乘法

% 那这个地方怎么去做呢? 怎么去考虑拓扑呢? 感觉难度比较大.

    coefs = [1;zeros(length(p_cumu),1)];
    for tempp = 1:length(p_cumu)
        % 这是一个主对角阵和一个次对角阵(下侧的)
        diag_matrix = diag((1-p_cumu(tempp)) * ones(1,length(p_cumu)+1),0) + diag( p_cumu(tempp) *  ones(1,length(p_cumu)),-1);
        coefs = diag_matrix * coefs; 
%         for i = tempp:-1:1
%             coefs(i+1) = coefs(i+1) * (1 - p_cumu(tempp)) + coefs(i) * p_cumu(tempp);
%         end 
%         coefs(1) = coefs(1) * (1 - p_cumu(tempp));
    end 
    
end
%概率性的,最终得到的是价格大于x的需求(供给)一共有多少个.
%我怎么得到对称的结果呢?似乎没有办法得到对称的结果.
%如何计算某个节点的剩余需求曲线(还得是概率性的),用来计算该节点是否能得到出清. 这个思路感觉根本就不可行. 因为维度太高,无法判断.
%如果我要内嵌一个优化问题,那我不如在原本的问题里去考虑.

% def f(ps):
%     coefs = [1]
%     for p in ps:
%         coefs.append(0)
%         for i in range(len(coefs) - 1, 0, -1):
%             coefs[i] = coefs[i] * (1 - p) + coefs[i - 1] * p
%         coefs[0] *= 1 - p
%     return coefs
   
       
    