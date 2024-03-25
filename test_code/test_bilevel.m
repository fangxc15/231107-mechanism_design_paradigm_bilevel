% 测试双层模型的求解效果

sdpvar x
sdpvar t
OO = t*t + 5*x;
CO = [t<=10];
OI = -x;
CI = [x<=t];
solvebilevel(CO,OO,CI,OI,x);
Value.t = value(t);
Value.x = value(x);

