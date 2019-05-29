function result=dwdt_w_movie(itwnt,w)

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
frmstep=10;
for i=itmin:frmstep:itmax
  ind=fsize*(i-1);

  plot(profs(1,ind+1:ind+fsize),profs(2,ind+1:ind+fsize),'-b');
  hold on
  %plot([w((i-1)/frmstep+1,2) w((i-1)/frmstep+1,2)],[-10 10],'r');
  %plot([w(i,2) w(i,2)],[-10 10],'r');
  plot(w(i,2),w(i,3),'ro');
  hold off
  ylabel('dw/dt');
  xlabel('w');
  axis([0 0.5 -20 20])

  f(i)=getframe;
%  mov=addframe(mov,f(i));
%  pause
end
%movie(f)
%mov=close(mov);
result=f;

