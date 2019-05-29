load p.d
[m n]=size(p);
x=min(p(:,1)):(max(p(:,1))-min(p(:,1)))/(10*(m-1)):max(p(:,1));
y=spline(p(:,1),p(:,2),x);
re=[x;y];
re=re';
save p.d2 re -ASCII
