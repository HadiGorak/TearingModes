function makeplots(nrm_in,colr)

eval(['!dwdt -db 0.1 -nrm ',num2str(nrm_in)]);
%-D_pol -0.24
load w.out
plot(w(:,1),w(:,2),'r');
hold on
!mv w.out w.out.0.1
bvsdbdt(1)=spline(w(:,2),w(:,5),[1.0]);
eval(['!dwdt -db 1 -nrm ',num2str(nrm_in)]);
load w.out
%plot(w(:,1),w(:,2),'g');
!mv w.out w.out.1.0
bvsdbdt(2)=spline(w(:,2),w(:,5),[1.0]);
eval(['!dwdt -db 10 -nrm ',num2str(nrm_in)]);
load w.out
%plot(w(:,1),w(:,2),'b');
!mv w.out w.out.10.0
bvsdbdt(3)=spline(w(:,2),w(:,5),[1.0]);
%hold off
%xlabel('time (s)');
%ylabel('W (cm)');

hold on
plot([0.1 1.0 10.0],bvsdbdt,'o-','Color',colr);
xlabel('d\beta_{N}/dt (s^{-1})');
ylabel('\beta(1cm)');
hold off
