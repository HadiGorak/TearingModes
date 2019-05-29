function [result]=makey(x)

[m n]=size(x);

for i=1:n
  result(i)=delta_prime(x(i));
end

