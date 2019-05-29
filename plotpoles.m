
num=3;

b_dbdt=zeros(num,1);
dwdt_b=zeros(num,1);
dwdt_w1=zeros(num,1);
w1=zeros(num,1);
dwdt_w2=zeros(num,1);
w2=zeros(num);
w_spec=0.5;
beta1=3.5;
beta2=3.6;
dbmin=2.0;
dbmax=15.0;
db=dbmin:(dbmax-dbmin)/(num-1):dbmax;
for i=1:num
  [b_dbdt,dwdt_b,dwdt_w1,w1,dwdt_w2,w2]=...
            plotonepole(i,'r',b_dbdt,dwdt_b,dwdt_w1,w1,dwdt_w2,w2,w_spec,...
            beta1,beta2,db(i));
  subplot(3,1,1);
  hold on
  subplot(3,1,2);
  hold on
  subplot(3,1,3);
  hold on
end

subplot(3,1,1);
ylabel('\beta_N');
hold off
axis([0 0.3 3.25 3.4]);
subplot(3,1,2);
ylabel('dp');
hold off
axis([0 0.3 -0.5 1.0]);
subplot(3,1,3);
ylabel('w');
xlabel('time (s)');
hold off
axis([0 0.3 0 5]);

figure
plot(db,b_dbdt);
xlabel('d\beta_{N}/dt (s^{-1})');
ylabel('\beta(1cm)');
print -deps bvsdbdt.eps

figure
plot(b_dbdt,dwdt_b);
xlabel('\beta(1cm)');
ylabel('dw/dt(1cm)');
print -deps dwdtvsbeta.eps

figure
plot(w1,dwdt_w1);
hold on
plot(w2,dwdt_w2,'r');
hold off
xlabel('w (cm)');
ylabel('dw/dt (cm/s)');
print -deps dwdtvsw.eps
w1(1)
