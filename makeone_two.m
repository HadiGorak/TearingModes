function makeone_two
%eval(['!dwdt -dtt 0.001 1.0 -db ',num2str(db),' -nrm ',num2str(nrm),...
%      ' -w0 ',num2str(w0)]);
!runner2
load w.out

subplot(5,1,1);
plot(w(:,1),w(:,6),'b');
ylabel('\beta_N');

subplot(5,1,2);
plot(w(:,1),w(:,7),'b');
ylabel('d\beta_N/dt');

subplot(5,1,3);
plot(w(:,1),w(:,5),'b');
ylabel('dp');

subplot(5,1,4);
plot(w(:,1),w(:,2),'b');
ylabel('w');
xlabel('time (s)');

subplot(5,1,5);
plot(w(:,1),w(:,4),'b');
ylabel('dwdt');
xlabel('time (s)');
%dwdt_w_movie(-1,w);

!runner3
load w.out

subplot(5,1,1);
hold on
plot(w(:,1),w(:,6),'r');
hold off

subplot(5,1,2);
hold on
plot(w(:,1),w(:,7),'r');
hold off

subplot(5,1,3);
hold on
plot(w(:,1),w(:,5),'r');
hold off

subplot(5,1,4);
hold on
plot(w(:,1),w(:,2),'r');
hold off

subplot(5,1,5);
hold on
plot(w(:,1),w(:,4),'r');
hold off
axis([0 1 -1 1]);
%dwdt_w_movie(-1,w);

