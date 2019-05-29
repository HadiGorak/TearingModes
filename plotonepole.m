function [b_dbdt,dwdt_b,dwdt_w1,w1,dwdt_w2,w2]=...
           plotonepole(index,col,b_dbdt,dwdt_b,dwdt_w1,w1,dwdt_w2,w2,w_spec, ...
                       beta1,beta2,db)

eval(['!dwdt -dtt 0.001 1.0 -db ',num2str(db),' -nrm 2.0 -w0 0.1']);
load w.out
subplot(3,1,1);
plot(w(:,1),w(:,5),col);
subplot(3,1,2);
plot(w(:,1),w(:,4),col);
subplot(3,1,3);
plot(w(:,1),w(:,2),col);

% find the index of the first change in w
[m n]=size(w);
k=1;
winit=w(1,2);
for i=1:m-1
  if w(i,2)==winit
    k=i;
  else
    break;
  end
end

% find the index of the first local maximum in w
% this may cause fail and if does, check pole peak condition
j=k;
for i=k+1:m-1
  if w(i,2)>w(j,2)
    j=i;
  else 
    break;
  end
end
b_dbdt(index)=spline(w(k:j,2),w(k:j,5),[w_spec]);
dwdt_b(index)=spline(w(k:j,2),w(k:j,3),[w_spec]);
dwdt_w1(index)=spline(w(k:j,5),w(k:j,3),[beta1]);
w1(index)=spline(w(k:j,5),w(k:j,2),[beta1]);
dwdt_w2(index)=spline(w(k:j,5),w(k:j,3),[beta2]);
w2(index)=spline(w(k:j,5),w(k:j,2),[beta2]);

