%function makeone(db,nrm,w0,col)
%eval(['!dwdt -dtt 0.001 1.0 -db ',num2str(db),' -nrm ',num2str(nrm),...
%      ' -w0 ',num2str(w0)]);

function runone(itwnt,db)
eval(['!runner4 ',num2str(db)]);

frmstep=5;

load w.out

fid=fopen('dwdt_w.out');
dumy=fscanf(fid,'%i',[1 1]);
mpsi=dumy(1);
profs=fscanf(fid,'%G %G',[2 inf]);
fclose(fid);

[m n]=size(profs);
fsize=mpsi;
itmax=n/fsize;

%mov = avifile('profs.avi')

if (itwnt==-1)
  itmin=1;
else
  itmin=itwnt;
  itmax=itwnt;
end

for i=itmin:frmstep:itmax
  ind=fsize*(i-1);

  subplot(4,1,1);
  plot(w(:,1),w(:,6),'b');
  hold on
  plot([w(i,1) w(i,1)],[1 2],'r');
  hold off
  ylabel('\beta_N');
  text(0.2,1.8,num2str(db));

  subplot(4,1,2);
  plot(w(:,1),w(:,5),'b');
  hold on
  plot([w(i,1) w(i,1)],[-5 5],'r');
  hold off
  ylabel('dp');

  subplot(4,1,3);
  plot(w(:,1),w(:,2),'b');
  hold on
  plot([w(i,1) w(i,1)],[0 0.3],'r');
  hold off
  ylabel('w');
  xlabel('time (s)');

  subplot(4,1,4);
  plot(profs(1,ind+1:ind+fsize),profs(2,ind+1:ind+fsize),'-b');
  hold on
  %plot([w((i-1)/frmstep+1,2) w((i-1)/frmstep+1,2)],[-10 10],'r');
  %plot([w(i,2) w(i,2)],[-10 10],'r');
  plot(w(i,2),w(i,4),'ro');
  hold off
  ylabel('dw/dt');
  xlabel('w');
  axis([0 0.5 -100 100])

  f(i)=getframe;
%  mov=addframe(mov,f(i));
%  pause
end
%movie(f)
%mov=close(mov);
result=f;

