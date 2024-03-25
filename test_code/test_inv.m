% 测试逆矩阵


Num.I = 2;
Num.P = 2;
Num.S = (Num.P + 1) ^Num.I;
% B = zeros(Num.S + Num.P * Num.N, Num.N * Num.S); %行是约束的个数

scenario_prob = ones(1,Num.S);
for i = 1:Num.I
    scenario_set(i,:) =  floor(mod(0:Num.S-1,(Num.P+1)^(Num.I- i +1))/((Num.P+1)^(Num.I-i)));
%     scenario_prob =  scenario_prob .* Para.Prob(i,scenario_set(i,:)+1);
end
B = repmat(diag(ones(Num.S,1)),1,Num.I);
for i = 1:Num.I
    for p = 0:Num.P
        temp = zeros(1,Num.I*Num.S);
        temp(find(scenario_set(i,:) == p) + (i-1) * Num.S) = 1;
        B = [B;temp];
    end
end
% B = B(1:size(B,1)-1,:);
rank(B) %都是不满秩的
invB = pinv(B); % 这个伪逆的求法是对的，这个可以得到伪逆前的若干系数
% invB = pinv(B'*B) * B'; % 逆矩阵是这样的
bleft = [zeros(Num.S,1); 0.8;0.7;0.6;-0.5;-0.7;-0.9];
Q = invB*bleft;
sum(Q .* Q) %这是求的二范数, 1.0467, 是最小的

%然后以后所有满足Bd = 0的新的d

rank(B)
rank(invB)
rank([B;Q']) % 相当于能和B垂直，那么一定能和Q垂直. 所以Q是一个二范数最小的量
%  这时候如果想让abs(b)下降
d = [1 0 -1 -1 0 1 0 0 0       -1 0 1 1 0 -1 0 0 0]';
Qnew = 0.01*d + Q;
max(abs(Qnew)) %从0.3333降低到了0.3233, 二范数会增大, 负无穷规范数会减小
B*Qnew
sum(Qnew .* Qnew)

%这个地方已经说明了,没有用,可以举出反例
%是不是要额外加一条约束? 也不能说明二范数最小啊. 可以尝试一下用二范数最小的效果和可行性。

%% 我再研究一下另外一种，写出Lagrange乘子的方式
height_b = size(B,1);
width_b  = size(B,2);

% 这个到底是什么矩阵?
% B矩阵是一个简单的限制矩阵

B_collect = [B zeros(height_b,height_b);  2*diag(ones(width_b,1)) B'];
bleft_collect = [bleft;zeros(width_b,1)];
resultb = pinv(B_collect) * bleft_collect; % 感觉这个地方会陷入一个满秩的情况. 所以要用pinv而不是inv
Q = resultb(1:width_b);
delta = resultb(width_b+1:width_b+height_b);%是这些约束的Lagrange乘子，无需求解可以直接得到这个问题. 这个Lagrange乘子是当右手项的约束变化的时候，我求得的最小二范数会怎么变,定义应该是这样.

% 如果改变的话，会导致最小的二范数损失多少。但是如果是带约束的优化问题的验证的话，那个没办法验证。所以只能这样。

% 假设说我知道想要改变多少并准备改变
%% 
% H = 2 * diag(ones(width_b,1));
R = pinv(B*B') * B;
R' * bleft;
% 也就是二次规划课件第5页的那个公式是对的
S = pinv(B*B');%直接写出封闭解
S*bleft %两倍的关系，可能有一点小的变化,一个是正负反掉了,一个是乘以2
%%
% invB = inv(B' * B) * B'; % 这个求伪逆的方法不一定对
% B * invB * B
% invB * B * invB

% 总的变量个数 2*3^2   = 18
% 总的约束个数 3^2+3*2 = 15; 如果说要是变量更多，那可能会 3^2 + 2*3*2 = 21
