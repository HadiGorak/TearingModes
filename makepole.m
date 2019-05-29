function result=makepole(norm,offset,maxy,betamin,betamax)
xx=1/1000:1/1000:1-1/1000;
yy=-pi*xx.*cot(pi*xx); 
xp=spline(yy,xx,[maxy]);
x=1/904:xp/904:xp;
y=spline(xx,yy,x);
x=betamin:(betamax-betamin)/902:betamax;
y=y+1.0;
y=y./(maxy+1);
y=(norm+offset).*y;
y=y-offset;
result(:,1)=x(:);
result(:,2)=y(:);
plot(x,y);       

