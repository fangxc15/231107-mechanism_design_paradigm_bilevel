function ConsOut = F_Complementarity_Yal(Cons, a, b, theta, Setting)
    %	本函数用于自动生成互补松弛条件
    %   本函数创建于20171120 16:57
    Cons = [Cons, a>=0];
    Cons = [Cons, b>=0];
%     if size(theta, 2) == 1
%         Cons = [Cons, implies(theta, a==0)];
%         Cons = [Cons, implies((1-theta), b==0)];
%     else
    if isfield(Setting,'Complementary_bigM')
        Cons = [Cons, a <= (1-theta) * Setting.Complementary_bigM];
        Cons = [Cons, b <= theta * Setting.Complementary_bigM];
    else
        Cons = [Cons, a <= (1-theta) * 1000];
        Cons = [Cons, b <= theta * 1000];
    end
    
%     for i = 1:size(theta,1)
%         for p = 1:size(theta,2)
% %             Cons = [Cons, implies(theta(:,p), a(:,p) == 0)];
% %             Cons = [Cons, implies(1-theta(:,p), b(:,p) == 0)];
% 
% %             Cons = [Cons, implies(theta(c,p), a(c,p) == 0)];
% %             Cons = [Cons, implies(1-theta(c,p), b(c,p) == 0)];
% 
%             Cons = [Cons, a(i,p) <= (1-theta(i,p)) * 1000];
%             Cons = [Cons, b(i,p) <= theta(i,p) * 1000];
%         end 
%     end
%     end
    ConsOut = Cons;
end