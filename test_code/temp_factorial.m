% 测试阶乘
%factorial
A = 3;
B = 2;
sum = 0;
for i = 1:B
    sum = sum + factorial(A) * factorial(B-1)/ factorial(A-i)/ factorial(i)/ factorial(B-i)/ factorial(i-1);
end 