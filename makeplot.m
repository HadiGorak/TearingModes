function makeplot

eval(['!dwdt -db 1 -dpolmod 2 -dpol -2e-5 -nrm 1.0 ']);
%-D_pol -0.24
load w.out
plot(w(:,1),w(:,2),'r');
hold on
!mv w.out w.out.0.1
eval(['!dwdt -db 10 -dpolmod 2 -dpol -2e-5 -nrm 1.0 ']);
load w.out
plot(w(:,1),w(:,2),'b');
!mv w.out w.out.10.0
bvsdbdt(3)=spline(w(:,2),w(:,5),[1.0]);
hold off
xlabel('time (s)');
ylabel('W (cm)');

