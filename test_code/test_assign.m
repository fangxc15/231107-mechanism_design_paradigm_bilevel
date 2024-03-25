% 测试assign

% t = sdpvar(1,1);
% x = [t,t];
% assign(t,4)
% value(x)

x= sdpvar(5,1);
xint = binvar(5,1);

Cons = [x>=[5;4;3;2;1].* xint,sum(xint)==1,[1,2,3,5,9]*x>=14];
% Cons = [x>=[5;4;3;2;1].* xint,sum(xint)==1];
% assign(x,[0,6,2,2,4]');
% assign(x(1),-2);
% assign(x(2),6);
assign(x(3),2);
assign(x(4),2);
assign(x(5),4);
optimize(Cons, sum([2;6;2;4;7].*x), sdpsettings('solver','gurobi','usex0',1));
% value(x)