
!runner2
load w.out

subplot(3,1,1);
plot(w(:,1),w(:,6),'r');
ylabel('\beta_N');

subplot(3,1,2);
plot(w(:,1),w(:,2),'r');
ylabel('w');
xlabel('time (s)');

subplot(3,1,3);
plot(w(:,1),w(:,7),'r');
ylabel('d\beta_N/dt');

!runner3
load w.out

subplot(3,1,1);
hold on
plot(w(:,1),w(:,6),'b');
hold off
%axis([0.1 0.6 1 2.5]);

subplot(3,1,2);
hold on
plot(w(:,1),w(:,2),'b');
hold off

subplot(3,1,3);
hold on
plot(w(:,1),w(:,7),'b');
hold off
axis([0 1 -15 5]);

%axis([0.1 0.6 0 0.3]);

