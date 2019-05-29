%function makeone(db,nrm,w0,col)
%eval(['!dwdt -dtt 0.001 1.0 -db ',num2str(db),' -nrm ',num2str(nrm),...
%      ' -w0 ',num2str(w0)]);

function makeone(col)
!runner4

load w.out
size(w)
subplot(4,1,1);
plot(w(:,1),w(:,6),col);
hold on
ylabel('\beta_N');

subplot(4,1,2);
plot(w(:,1),w(:,5),col);
ylabel('dp');

subplot(4,1,3);
plot(w(:,1),w(:,2),col);
ylabel('w');
xlabel('time (s)');

subplot(4,1,4);
%plot(w(:,1),w(:,4),col);
%ylabel('dwdt');
%xlabel('time (s)');
dwdt_w_movie(-1,w);

subplot(4,1,1);
hold off
subplot(4,1,2);
hold off
subplot(4,1,3);
hold off
