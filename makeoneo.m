
!dwdt -dtt 0.0001 0.5 -db 2.0 -nrm 1.0
load w.out
subplot(3,1,1);
plot(w(:,1),w(:,5));
subplot(3,1,2);
plot(w(:,1),w(:,4));
subplot(3,1,3);
plot(w(:,1),w(:,2));
hold on
plot(w(:,1),w(:,2),'o');
hold off

